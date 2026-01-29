import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Path to the tracked CSV file (do not change; this points to CardilloMPI/vtk_out)
csv_path = '../../../../latex/plots/data/woodpecker_tracked.csv'

def main():
    # Read CSV
    df = pd.read_csv(csv_path)

    # Ensure 't' exists
    time = df['t'] if 't' in df else df.iloc[:,0]

    # Print a few lines for quick inspection
    print('CSV head:')
    print(df.head())

    # Filter rows for ring and body
    df_ring = df[df['name'] == 'woodpecker_ring']
    df_body = df[df['name'] == 'woodpecker_body']

    # Create 2 stacked subplots (body on top, ring on bottom)
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    fig.suptitle('Woodpecker: Body (top) and Ring (bottom)')

    # Body subplot: primary y for pz/vz, secondary y for wz and phi_z (euler_z)
    ax = axes[0]
    if df_body.empty:
        ax.text(0.5, 0.5, 'No data for woodpecker_body', ha='center')
    else:
        if 'pz' in df_body:
            ax.plot(df_body['t'], df_body['pz'], label='pz')
        if 'vz' in df_body:
            ax.plot(df_body['t'], df_body['vz'], label='vz')
        ax.set_ylabel('pz / vz')
        ax.set_title('Body: pz, vz (left) and wz, phi_z (right)')
        ax.grid(True)
        ax.legend(loc='upper left')

        ax2 = ax.twinx()
        if 'wz' in df_body:
            ax2.plot(df_body['t'], df_body['wz'], color='C3', linestyle='--', label='wz')
        if 'euler_z' in df_body:
            ax2.plot(df_body['t'], df_body['euler_z'], color='C4', linestyle=':', label='phi_z')
        ax2.set_ylabel('wz / phi_z')
        # combine legends
        lines, labels = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2, loc='upper right')

    # Ring subplot: primary y for pz/vz, secondary y for wy and phi_y (euler_y)
    ax = axes[1]
    if df_ring.empty:
        ax.text(0.5, 0.5, 'No data for woodpecker_ring', ha='center')
    else:
        if 'pz' in df_ring:
            ax.plot(df_ring['t'], df_ring['pz'], label='pz')
        if 'vz' in df_ring:
            ax.plot(df_ring['t'], df_ring['vz'], label='vz')
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('pz / vz')
        ax.set_title('Ring: pz, vz (left) and wy, phi_y (right)')
        ax.grid(True)
        ax.legend(loc='upper left')

        ax2 = ax.twinx()
        if 'wy' in df_ring:
            ax2.plot(df_ring['t'], df_ring['wy'], color='C3', linestyle='--', label='wy')
        if 'euler_y' in df_ring:
            ax2.plot(df_ring['t'], df_ring['euler_y'], color='C4', linestyle=':', label='phi_y')
        ax2.set_ylabel('wy / phi_y')
        lines, labels = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2, loc='upper right')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

if __name__ == '__main__':
    main()
