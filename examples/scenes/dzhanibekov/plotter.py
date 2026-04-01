import os

csv_folder = "../../../vtk_out"
# trim accidental whitespace
csv_filename = "dzhanibekov_tracked.csv"

#t,name,px,py,pz,vx,vy,vz,wx,wy,wz,euler_x,euler_y,euler_z
#plot tracked euler angles over time

import matplotlib.pyplot as plt
import pandas as pd



csv_path = os.path.normpath(os.path.join(os.path.dirname(__file__), csv_folder, csv_filename))
if not os.path.exists(csv_path):
	raise SystemExit(f"Tracked CSV not found: {csv_path}\nRun the simulation (main) to produce vtk_out, or point this script to the correct file.")

df = pd.read_csv(csv_path)
plt.figure(figsize=(10,6))
plt.plot(df['t'], df['wx'], color='red', label='wx')
plt.plot(df['t'], df['wy'], color='blue', label='wy')
plt.plot(df['t'], df['wz'], color='black', label='wz')
plt.xlabel('Time (s)')

plt.ylabel('Angular velocity (rad/s)')
plt.title('Tracked Angular Velocities Over Time')
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()