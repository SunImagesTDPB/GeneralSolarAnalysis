import numpy as np
from scipy.interpolate import interp1d
import os
import sunpy.map
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
import matplotlib.pyplot as plt
from astropy.visualization import ImageNormalize, AsinhStretch

lines = [94, 131, 171, 193, 211, 304, 335]


def line_to_pixels_with_values(image, grid_x, grid_y, control_x, control_y):

    dx = np.diff(np.round(image.shape[1]*(control_x-grid_x[0])/(grid_x[-1]-grid_x[0])))[0]
    dy = np.diff(np.round(image.shape[0]*(control_y-grid_y[0])/(grid_y[-1]-grid_y[0])))[0]
    num_points = np.round(np.sqrt(dx**2+dy**2)).astype(int)
    t = np.linspace(0, 1, len(control_x))
    t_new = np.linspace(0, 1, num_points)
    lx = interp1d(t, control_x)
    ly = interp1d(t, control_y)
    x = lx(t_new)
    y = ly(t_new)
    x_pixels = np.round(image.shape[1]*(x-grid_x[0])/(grid_x[-1]-grid_x[0])).astype(int)
    y_pixels = np.round(image.shape[0]*(y-grid_y[0])/(grid_y[-1]-grid_y[0])).astype(int)
    
    pixel_values = []
    for i in range(num_points):
        if 0 <= x_pixels[i] < image.shape[1] and 0 <= y_pixels[i] < image.shape[0]:
            pixel_values.append(image[y_pixels[i], x_pixels[i]])
        else:
            pixel_values.append(None)
    
    return pixel_values


def gen_slices(x_loop, y_loop, width):
    ts = np.linspace(0,1,len(x_loop))
    xs = CubicSpline(ts, x_loop)
    ys = CubicSpline(ts, y_loop)
    tanvecs = np.array([np.array([xs(ts[i],1), ys(ts[i],1)])/np.linalg.norm(np.array([xs(ts[i],1), ys(ts[i],1)])) for i in range(len(x_loop))])
    normvecs = np.array([tanvecs[:,1], -tanvecs[:,0]]).T
    slicex = np.array([x_loop-0.5*width*normvecs[:,0], x_loop+0.5*width*normvecs[:,0]]).T
    slicey = np.array([y_loop-0.5*width*normvecs[:,1], y_loop+0.5*width*normvecs[:,1]]).T
    return slicex, slicey


class TDPlotter:
    def __init__(self, dir, channel=131):
        images = os.listdir(dir)
        self.aia_maps = {}
        self.channel = channel
        for l in lines:
            tmp = []
            for filename in images:
                if filename.find(f'{l}A') != -1:
                    tmp.append(sunpy.map.Map(f'{dir}/{filename}'))
            self.aia_maps[l] = sunpy.map.Map(tmp, sequence=True)
        self.main_map = self.aia_maps[self.channel][0]
        self.x_grid = np.linspace(self.main_map.bottom_left_coord.Tx.value, self.main_map.top_right_coord.Tx.value, self.main_map.data.shape[1])
        self.y_grid = np.linspace(self.main_map.bottom_left_coord.Ty.value, self.main_map.top_right_coord.Ty.value, self.main_map.data.shape[0])
        self.time_array = np.array([m.date.tai_seconds for m in self.aia_maps[self.channel]])-self.main_map.date.tai_seconds

    def get_data(self, x_loop, y_loop, width):
        self.dist = np.zeros(len(x_loop))
        self.dist[1:] = np.cumsum(np.sqrt(np.diff(x_loop)**2+np.diff(y_loop)**2))
        self.width = width
        slicex, slicey = gen_slices(x_loop, y_loop, width)
        self.intensity_stack = {}
        for l in lines:
            self.intensity_stack[str(l)] = []
            for i in range(len(self.aia_maps[l])):
                prevalues = [line_to_pixels_with_values(self.aia_maps[l][i].data, self.x_grid, self.y_grid, slicex[j], slicey[j]) for j in range(len(x_loop))]
                self.intensity_stack[str(l)].append([np.mean(v) for v in prevalues])
            self.intensity_stack[str(l)] = np.array(self.intensity_stack[str(l)])
    
    def plot_td(self, channel=131):
        fig, ax = plt.subplots(figsize=(15, 5))
        plt.pcolormesh(self.time_array, self.dist, self.intensity_stack[str(channel)].T,
                    norm=ImageNormalize(stretch=AsinhStretch()), cmap=f'sdoaia{channel}')
        ax.set_ylabel("Distance along slit (arcsec)")
        ax.set_xlabel("Time (s)")
        plt.show()
    
    def save_npz(self, filename):
        self.intensity_stack['dist'] = self.dist * 0.72527094 # перевод длины петли в мегаметры
        self.intensity_stack['w'] = self.width * 0.72527094 # перевод ширины петли в мегаметры
        np.savez(f'{filename}.npz', **self.intensity_stack)
    
    def save_txt(self, filename):
        np.savetxt(f'{filename}_x_grid.txt', self.dist)
        np.savetxt(f'{filename}_t_grid.txt', self.time_array)
        for l in lines:
            np.savetxt(f'{filename}_{l}_line.txt', self.intensity_stack[str(l)])
