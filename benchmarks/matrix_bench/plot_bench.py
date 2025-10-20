#!/usr/bin/env python3
import os
os.environ.setdefault('MPLBACKEND', 'Agg')  # force non-interactive backend
import subprocess
import re
import argparse
import numpy as np
import matplotlib
matplotlib.use(os.environ.get('MPLBACKEND', 'Agg'))
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BIN_CANDIDATES = [
    os.path.normpath(os.path.join(SCRIPT_DIR, '../../build/bin/matrix_bench')),
    os.path.normpath(os.path.join(SCRIPT_DIR, '../..', 'build', 'bin', 'matrix_bench')),
]
BIN = None
for c in BIN_CANDIDATES:
    if os.path.exists(c):
        BIN = c
        break
if BIN is None:
    raise SystemExit("matrix_bench binary not found. Build the project first.")

def run_once(impl, nb, nc, iters=5, rigid_frac=0.5):
    cmd = [BIN, '--impl', impl, '--nb', str(nb), '--nc', str(nc), '--iters', str(iters), '--rigid-frac', str(rigid_frac)]
    out = subprocess.check_output(cmd, text=True)
    # parse timings from stdout
    tm = {}
    for line in out.splitlines():
        m = re.search(r'A\*x:\s+([0-9.]+) ms', line)
        if m: tm['mul'] = float(m.group(1))
        m = re.search(r'rowDot sum:\s+([0-9.]+) ms', line)
        if m: tm['row'] = float(m.group(1))
        m = re.search(r'At\*w \(acc\):\s+([0-9.]+) ms', line)
        if m: tm['atw'] = float(m.group(1))
        m = re.search(r'At\*A\*x:\s+([0-9.]+) ms', line)
        if m: tm['normal'] = float(m.group(1))
        m = re.search(r'sum\(y\)=([\-0-9.eE]+)', line)
        if m: tm['sum_y'] = float(m.group(1))
        m = re.search(r'sum\(z\)=([\-0-9.eE]+)', line)
        if m: tm['sum_z'] = float(m.group(1))
        m = re.search(r'sum\(g\)=([\-0-9.eE]+)', line)
        if m: tm['sum_g'] = float(m.group(1))
        if 'Verify A*x equal:' in line:
            tm['verify'] = line.strip()
    return tm

def _close(a, b, tol=1e-6):
    if a is None or b is None:
        return False
    da = abs(a - b)
    m = max(1.0, abs(a) + abs(b))
    return da <= 1e-6 or da/m <= tol

def validate_equivalent(nb, nc, impls, rigid_frac=0.5, tol=1e-6):
    # Run each impl once and compare sums of y and z to the first impl
    order = list(impls)
    results = {}
    for impl in order:
        tm = run_once(impl, nb, nc, iters=1, rigid_frac=rigid_frac)
        results[impl] = tm
    ref = results[order[0]]
    ok = True
    for impl, tm in results.items():
        if (not _close(tm.get('sum_y'), ref.get('sum_y'), tol)
            or not _close(tm.get('sum_z'), ref.get('sum_z'), tol)
            or not _close(tm.get('sum_g'), ref.get('sum_g'), tol)):
            print(f"Validation mismatch nb={nb} impl={impl}: sums differ\n  ref: y={ref.get('sum_y')} g={ref.get('sum_g')} z={ref.get('sum_z')}\n  got: y={tm.get('sum_y')} g={tm.get('sum_g')} z={tm.get('sum_z')}")
            ok = False
    if ok:
        print(f"Validation OK at nb={nb} across {', '.join(order)}")
    return ok

def sweep_sizes(impls=('block','sparse','dense','blas','cublas'), nb_list=(500,1000,2000), density=2.0, iters=5, rigid_frac=0.5):
    results = {impl: {'nb': [], 'nc': [], 'mul': [], 'row': [], 'atw': [], 'normal': []} for impl in impls}
    for nb in nb_list:
        nc = int(density*nb)
        # Validate numerical equivalence on a subset including GPU backends
        try:
            subset = [i for i in impls if i in ('block','sparse','cublas','cusparse')]
            # if subset:
            #     validate_equivalent(nb, nc, subset, rigid_frac=rigid_frac)
        except Exception as e:
            print(f"Validation step failed at nb={nb}: {e}")
        for impl in impls:
            if (impl == 'dense' and nb > 2000) or (impl == 'blas' and nb > 2000) or (impl == 'cublas' and nb > 5000):
                # skip dense for very large sizes to keep runtime sane
                for k in ['mul','row','atw','normal']:
                    results[impl][k].append(np.nan)
                results[impl]['nb'].append(nb)
                results[impl]['nc'].append(nc)
                continue
            try:
                tm = run_once(impl, nb, nc, iters=iters, rigid_frac=rigid_frac)
            except subprocess.CalledProcessError as e:
                print(f"Skipping impl '{impl}' at nb={nb}: {e}")
                for k in ['mul','row','atw','normal']:
                    results[impl][k].append(np.nan)
                results[impl]['nb'].append(nb)
                results[impl]['nc'].append(nc)
                continue
            for k in ['mul','row','atw','normal']:
                results[impl][k].append(tm.get(k, np.nan))
            results[impl]['nb'].append(nb)
            results[impl]['nc'].append(nc)
            print(f"{impl} nb={nb} nc={nc} -> {tm}")
    return results

def plot_results(results, title='Matrix ops benchmark (nc ≈ 2·nb)', save='bench_plot.png'):
    fig, axes = plt.subplots(2, 2, figsize=(12,8), sharex=True, sharey=True)
    ops = ['mul','row','atw','normal']
    op_titles = {'mul':'A*x', 'row':'row dot sum', 'atw':'A^T*w (acc)', 'normal':'A^T*A*x'}
    styles = {
        'block': {'color':'tab:blue', 'ls':'-', 'marker':'o', 'label':'BlockArray'},
        'sparse': {'color':'tab:green', 'ls':'--', 'marker':'s', 'label':'EigenSparse'},
        'dense': {'color':'tab:red', 'ls':':', 'marker':'D', 'label':'EigenDense'},
        'blas': {'color':'tab:orange', 'ls':'-.', 'marker':'^', 'label':'BLAS'},
        'cublas': {'color':'tab:purple', 'ls':'-', 'marker':'v', 'label':'cuBLAS'},
        'cusparse': {'color':'tab:brown', 'ls':'--', 'marker':'x', 'label':'cuSPARSE'},
    }
    for ax, op in zip(axes.flatten(), ops):
        for impl, data in results.items():
            x = data['nb']; y = data[op]
            if len(x)==0: continue
            s = styles.get(impl, {'color':'k','ls':'-','marker':'.','label':impl})
            ax.plot(x, y, linestyle=s['ls'], color=s['color'], marker=s['marker'], label=s['label'])
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_title(op_titles[op])
        ax.grid(True, which='both', ls='--', alpha=0.3)
    handles, labels = axes[0,0].get_legend_handles_labels()
    bylabel = dict(zip(labels, handles))
    fig.legend(bylabel.values(), bylabel.keys(), loc='upper center', ncol=5, fontsize=9)
    fig.suptitle(title)
    fig.text(0.5, 0.04, 'num bodies (nb)', ha='center')
    fig.text(0.04, 0.5, 'time (ms)', va='center', rotation='vertical')
    plt.tight_layout(rect=[0.02, 0.05, 0.98, 0.90])
    plt.savefig(save, dpi=200)
    plt.close(fig)
    print(f"Saved plot to {save}")

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--min-nb', type=int, default=500)
    ap.add_argument('--max-nb', type=int, default=2000)
    ap.add_argument('--step', type=int, default=500)
    ap.add_argument('--density', type=float, default=2.0)
    ap.add_argument('--iters', type=int, default=5)
    ap.add_argument('--rigid-frac', type=float, default=0.5)
    ap.add_argument('--save', type=str, default='bench_plot.png')
    ap.add_argument('--impls', type=str, default='block,sparse,dense,blas,cublas,cusparse')
    args = ap.parse_args()

    nb_list = tuple(range(args.min_nb, args.max_nb+1, args.step))
    impls = tuple([s.strip() for s in args.impls.split(',') if s.strip()])
    results = sweep_sizes(impls=impls, nb_list=nb_list, density=args.density, iters=args.iters, rigid_frac=args.rigid_frac)
    plot_results(results, save=args.save)