#!/usr/bin/env python3

import argparse
import math
import os
import shutil
import sys
import xml.etree.ElementTree as ET
import numpy as np

try:
    import vtk
    from vtk.util import numpy_support as vnp
except ImportError:
    print(
        "This script requires the 'vtk' python package.\nInstall with:  pip install vtk numpy",
        file=sys.stderr,
    )
    raise
XML_READER_WRITER = {
    ".vtu": (vtk.vtkXMLUnstructuredGridReader, vtk.vtkXMLUnstructuredGridWriter),
    ".vtp": (vtk.vtkXMLPolyDataReader, vtk.vtkXMLPolyDataWriter),
    ".vts": (vtk.vtkXMLStructuredGridReader, vtk.vtkXMLStructuredGridWriter),
    ".vtr": (vtk.vtkXMLRectilinearGridReader, vtk.vtkXMLRectilinearGridWriter),
    ".vti": (vtk.vtkXMLImageDataReader, vtk.vtkXMLImageDataWriter),
}


def get_reader_writer(ext):
    ext = ext.lower()
    if ext == ".vtk":
        return (vtk.vtkDataSetReader, vtk.vtkDataSetWriter)
    if ext in XML_READER_WRITER:
        return XML_READER_WRITER[ext]
    raise ValueError(
        f"Unsupported VTK data file extension: '{ext}'. Supported: .vtk, {', '.join(XML_READER_WRITER)}"
    )


def read_dataset(path):
    ext = os.path.splitext(path)[1]
    reader_cls, _ = get_reader_writer(ext)
    reader = reader_cls()
    reader.SetFileName(path)
    reader.Update()
    return (reader.GetOutput(), ext)


def write_dataset(dataset, path, ext, binary=True):
    _, writer_cls = get_reader_writer(ext)
    writer = writer_cls()
    writer.SetFileName(path)
    writer.SetInputData(dataset)
    if binary and hasattr(writer, "SetDataModeToBinary"):
        writer.SetDataModeToBinary()
    elif binary and hasattr(writer, "SetFileTypeToBinary"):
        try:
            writer.SetFileTypeToBinary()
        except Exception:
            pass
    writer.Write()


def get_point_array(dataset, name):
    arr = dataset.GetPointData().GetArray(name)
    if arr is None:
        raise ValueError(
            f"Point-data array '{name}' not found. Available arrays: {[dataset.GetPointData().GetArrayName(i) for i in range(dataset.GetPointData().GetNumberOfArrays())]}"
        )
    return vnp.vtk_to_numpy(arr)


def parse_pvd(pvd_path):
    tree = ET.parse(pvd_path)
    root = tree.getroot()
    collection = root.find("Collection")
    if collection is None:
        raise ValueError(
            f"'{pvd_path}' does not look like a valid ParaView .pvd Collection file (no <Collection> element)."
        )
    entries = []
    for ds in collection.findall("DataSet"):
        ts = ds.get("timestep")
        entries.append(
            {"timestep": float(ts) if ts is not None else None, "file": ds.get("file")}
        )
    return entries


def chronological_order(entries):
    return sorted(
        range(len(entries)),
        key=lambda i: (
            entries[i]["timestep"] if entries[i]["timestep"] is not None else i
        ),
    )


def parse_aspect_ratio(s):
    if ":" in s:
        w, h = s.split(":")
        return float(w) / float(h)
    return float(s)


def rotation_matrix_2d(deg):
    th = math.radians(deg)
    c, s = (math.cos(th), math.sin(th))
    return np.array([[c, -s], [s, c]])


def extract_plane_columns(points, plane):
    idx = {"xy": (0, 1), "xz": (0, 2), "yz": (1, 2)}[plane]
    return points[:, idx]


def compute_tex_coords(
    points_2d,
    rotation_deg,
    center_override,
    scale_mult,
    aspect_ratio,
    fit_margin,
    flip_u,
    flip_v,
):
    pts = points_2d.astype(np.float64).copy()
    if rotation_deg % 360.0 != 0.0:
        pts = pts @ rotation_matrix_2d(rotation_deg).T
    if center_override is not None:
        cx, cy = center_override
    else:
        cx = (pts[:, 0].min() + pts[:, 0].max()) / 2.0
        cy = (pts[:, 1].min() + pts[:, 1].max()) / 2.0
    dx = pts[:, 0] - cx
    dy = pts[:, 1] - cy
    bw = dx.max() - dx.min()
    bh = dy.max() - dy.min()
    bw = bw if bw > 1e-12 else 1.0
    bh = bh if bh > 1e-12 else 1.0
    span_x = bw
    span_y = bh * aspect_ratio
    base = max(span_x, span_y)
    s = fit_margin / base if base > 0 else 1.0
    s *= scale_mult
    u = 0.5 + dx * s
    v = 0.5 + dy * s * aspect_ratio
    if flip_u:
        u = 1.0 - u
    if flip_v:
        v = 1.0 - v
    u = np.clip(u, 0.0, 1.0)
    v = np.clip(v, 0.0, 1.0)
    uv = np.stack([u, v], axis=1).astype(np.float32)
    return (uv, (cx, cy), s)


def build_entity_uv_map(entity_ids, uv):
    mapping = {}
    counts = {}
    for eid in np.unique(entity_ids):
        idxs = np.nonzero(entity_ids == eid)[0]
        mapping[eid] = uv[idxs]
        counts[eid] = len(idxs)
    sizes = set(counts.values())
    if len(sizes) > 1:
        print(
            f"    note: reference frame entities have varying vertex counts {sorted(sizes)}; expected a single fixed count."
        )
    return (mapping, counts)


def assign_uv_by_entity(entity_ids, mapping):
    n = len(entity_ids)
    uv_out = np.zeros((n, 2), dtype=np.float32)
    missing = 0
    mismatched = 0
    for eid in np.unique(entity_ids):
        idxs = np.nonzero(entity_ids == eid)[0]
        ref_uv = mapping.get(eid)
        if ref_uv is None:
            missing += 1
            continue
        k = min(len(idxs), len(ref_uv))
        if len(idxs) != len(ref_uv):
            mismatched += 1
        uv_out[idxs[:k]] = ref_uv[:k]
    return (uv_out, missing, mismatched)


def set_tex_coords(dataset, uv, array_name):
    n = dataset.GetNumberOfPoints()
    if uv.shape[0] != n:
        raise ValueError(
            f"Point count mismatch: dataset has {n} points but the tex-coord array being written has {uv.shape[0]}."
        )
    vtk_arr = vnp.numpy_to_vtk(np.ascontiguousarray(uv, dtype=np.float32), deep=True)
    vtk_arr.SetName(array_name)
    pd = dataset.GetPointData()
    if pd.HasArray(array_name):
        pd.RemoveArray(array_name)
    pd.AddArray(vtk_arr)
    pd.SetTCoords(vtk_arr)


def main():
    ap = argparse.ArgumentParser(
        description="Bake texture coordinates from one reference frame's vertex positions into a ParaView .pvd time series (for a 'rocks fall into a picture' effect).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("pvd", help="Path to the input .pvd file")
    ap.add_argument(
        "--mode",
        choices=["single", "all"],
        default="all",
        help="'single': tag ONLY the reference frame's file (fast -- no other files are read or written). 'all' (default): also propagate the same per-entity tex coords, matched by entity_id, into every other frame's file.",
    )
    ap.add_argument(
        "--frame-index",
        type=int,
        default=-1,
        help="Chronological index of the reference frame (Python-style, -1 = last). Default: -1.",
    )
    ap.add_argument(
        "--frame-time",
        type=float,
        default=None,
        help="Pick the reference frame by timestep value (closest match) instead of --frame-index.",
    )
    ap.add_argument(
        "--list-frames",
        action="store_true",
        help="Print available frames (index, timestep, file) and exit.",
    )
    ap.add_argument(
        "--entity-id-array",
        type=str,
        default="entity_id",
        help="Name of the point-data array giving each vertex's (unique per entity) entity id. Default: 'entity_id'.",
    )
    ap.add_argument(
        "--plane",
        choices=["xy", "xz", "yz"],
        default="xy",
        help="Which two world axes to project onto the image. Default: xy.",
    )
    ap.add_argument(
        "--rotation",
        type=float,
        default=0.0,
        help="Rotate the projected points by this many degrees about the center before mapping to UV (e.g. 90, 180, 270).",
    )
    ap.add_argument(
        "--center",
        type=float,
        nargs=2,
        default=None,
        metavar=("CX", "CY"),
        help="Override the auto-computed center (world units, measured AFTER rotation). Default: bounding-box center of the reference frame's points.",
    )
    ap.add_argument(
        "--scale",
        type=float,
        default=1.0,
        help="Extra zoom multiplier on top of the automatic 'fit the whole shape into frame' scale. >1 zooms in, <1 zooms out. Default: 1.0.",
    )
    ap.add_argument(
        "--aspect-ratio",
        type=str,
        default="1:1",
        help="Target image aspect ratio, width:height (e.g. '16:9') or a decimal (e.g. '1.778'). Default: '1:1' (square).",
    )
    ap.add_argument(
        "--fit-margin",
        type=float,
        default=0.95,
        help="Fraction of the [0,1] UV range the auto-fit scale should occupy, i.e. padding around the shape. Default: 0.95.",
    )
    ap.add_argument("--flip-u", action="store_true", help="Mirror U (horizontal flip).")
    ap.add_argument(
        "--flip-v",
        dest="flip_v",
        action="store_true",
        default=True,
        help="Flip V (default ON: image row 0 is usually the top, while world Y usually points up).",
    )
    ap.add_argument(
        "--no-flip-v",
        dest="flip_v",
        action="store_false",
        help="Disable the default V flip.",
    )
    ap.add_argument(
        "--array-name",
        type=str,
        default="TextureCoordinates",
        help="Name of the point-data array to write (also set as the active TCoords). Default: 'TextureCoordinates'.",
    )
    ap.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Write output file(s) (+ new .pvd) here instead of overwriting in place. In --mode all, all frames are copied here. In --mode single, only the reference frame is written here (use this if you don't want to touch the original file at all).",
    )
    ap.add_argument(
        "--backup",
        action="store_true",
        help="When overwriting in place, save each modified file as <file>.bak first.",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Compute and print bounding box / center / scale info, but write nothing.",
    )
    args = ap.parse_args()
    pvd_path = os.path.abspath(args.pvd)
    pvd_dir = os.path.dirname(pvd_path)
    entries = parse_pvd(pvd_path)
    if not entries:
        sys.exit("No <DataSet> entries found in the .pvd file.")
    order = chronological_order(entries)
    if args.list_frames:
        print(f"{'idx':>4}  {'timestep':>12}  file")
        for pos, orig_i in enumerate(order):
            e = entries[orig_i]
            ts = e["timestep"] if e["timestep"] is not None else pos
            print(f"{pos:>4}  {ts:>12}  {e['file']}")
        return
    if args.frame_time is not None:
        best_pos = min(
            range(len(order)),
            key=lambda p: abs(
                (
                    entries[order[p]]["timestep"]
                    if entries[order[p]]["timestep"] is not None
                    else p
                )
                - args.frame_time
            ),
        )
        ref_orig_i = order[best_pos]
    else:
        ref_orig_i = order[args.frame_index]
    ref_entry = entries[ref_orig_i]
    ref_path = os.path.join(pvd_dir, ref_entry["file"])
    print(f"Reference frame: {ref_entry['file']} (timestep={ref_entry['timestep']})")
    ref_dataset, ref_ext = read_dataset(ref_path)
    n_points = ref_dataset.GetNumberOfPoints()
    if n_points == 0:
        sys.exit(f"Reference frame '{ref_path}' has zero points.")
    ref_entity_ids = get_point_array(ref_dataset, args.entity_id_array)
    n_entities = len(np.unique(ref_entity_ids))
    points_np = vnp.vtk_to_numpy(ref_dataset.GetPoints().GetData())
    plane_pts = extract_plane_columns(points_np, args.plane)
    aspect_ratio = parse_aspect_ratio(args.aspect_ratio)
    uv, center_used, scale_used = compute_tex_coords(
        plane_pts,
        rotation_deg=args.rotation,
        center_override=args.center,
        scale_mult=args.scale,
        aspect_ratio=aspect_ratio,
        fit_margin=args.fit_margin,
        flip_u=args.flip_u,
        flip_v=args.flip_v,
    )
    print(
        f"Points: {n_points}   Entities: {n_entities}   Center used: {center_used}   Effective scale: {scale_used:.6g}   Aspect ratio (W/H): {aspect_ratio:.4g}"
    )
    print(
        f"U range: [{uv[:, 0].min():.4f}, {uv[:, 0].max():.4f}]   V range: [{uv[:, 1].min():.4f}, {uv[:, 1].max():.4f}]"
    )
    if args.dry_run:
        print("Dry run: no files were written.")
        return
    out_dir = os.path.abspath(args.output_dir) if args.output_dir else None
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    def out_path_for(rel_file):
        if out_dir:
            p = os.path.join(out_dir, rel_file)
            os.makedirs(os.path.dirname(p) or ".", exist_ok=True)
            return p
        return os.path.join(pvd_dir, rel_file)

    if args.mode == "single":
        dst_path = out_path_for(ref_entry["file"])
        set_tex_coords(ref_dataset, uv, args.array_name)
        if not out_dir and args.backup and os.path.exists(ref_path):
            shutil.copy2(ref_path, ref_path + ".bak")
        write_dataset(ref_dataset, dst_path, ref_ext)
        print(
            f"Wrote texture coordinates into 1 file ({dst_path}). No other frames were touched."
        )
        if out_dir:
            new_pvd_path = os.path.join(out_dir, os.path.basename(pvd_path))
            shutil.copy2(pvd_path, new_pvd_path)
            print(
                f"NOTE: only '{ref_entry['file']}' was written to '{out_dir}'. Copied the original .pvd there too, but the other frame files referenced by it were NOT copied (mode=single skips them to save time)."
            )
        return
    mapping, ref_counts = build_entity_uv_map(ref_entity_ids, uv)

    def process_file(rel_file):
        src_path = os.path.join(pvd_dir, rel_file)
        dst_path = out_path_for(rel_file)
        if src_path == ref_path:
            dataset, ext = (ref_dataset, ref_ext)
            set_tex_coords(dataset, uv, args.array_name)
        else:
            dataset, ext = read_dataset(src_path)
            entity_ids = get_point_array(dataset, args.entity_id_array)
            uv_out, missing, mismatched = assign_uv_by_entity(entity_ids, mapping)
            if missing:
                print(
                    f"  note: '{rel_file}': {missing} entity id(s) not found in the reference frame, left at UV (0, 0)."
                )
            if mismatched:
                print(
                    f"  note: '{rel_file}': {mismatched} entit(y/ies) had a different vertex count than in the reference frame; truncated/padded."
                )
            set_tex_coords(dataset, uv_out, args.array_name)
        if not out_dir and args.backup and os.path.exists(src_path):
            shutil.copy2(src_path, src_path + ".bak")
        write_dataset(dataset, dst_path, ext)
        return True

    count = 0
    for e in entries:
        if process_file(e["file"]):
            count += 1
    print(
        f"Wrote texture coordinates into {count}/{len(entries)} files{(' (output dir: ' + out_dir + ')' if out_dir else ' (in place)')}."
    )
    if out_dir:
        new_pvd_path = os.path.join(out_dir, os.path.basename(pvd_path))
        shutil.copy2(pvd_path, new_pvd_path)
        print(f"Copied .pvd to: {new_pvd_path}")


if __name__ == "__main__":
    main()
