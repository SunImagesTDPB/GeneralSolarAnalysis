from IPython.display import display
import ipywidgets as widgets
from ipywidgets import Output, Button, VBox

import numpy as np
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt

import sunpy.map
import sunpy.data.sample
from astropy.io import fits
from astropy.visualization import AsymmetricPercentileInterval, ImageNormalize, LogStretch, AsinhStretch

from sunpy.net import Fido, attrs as a
from astropy import units as u

import os
import time
import tkinter as tk

class InteractiveLoopTracer:
    def __init__(self, image_data, x_grid, y_grid):
        self.image_data = image_data
        self.x_grid = x_grid
        self.y_grid = y_grid
        self.fig, self.ax = plt.subplots(figsize=(4, 6))
        # Добавил нормализацию картинки, чтобы можно было разобрать петлю в нужном канале
        self.img = self.ax.pcolormesh(self.x_grid, self.y_grid, self.image_data, #origin='lower',
                                  norm=ImageNormalize(stretch=AsinhStretch()), cmap=f'sdoaia131')
        # self.img.set_clim(vmin=np.percentile(image_data, 1),
        #                   vmax=np.percentile(image_data, 99))

        self.points = []
        self.point_artists = []
        self.line_artist = None

        self.output = widgets.Output()
        self.coordinates_output = widgets.Output()

        self.clear_button = widgets.Button(description="Очистить точки")
        self.interpolate_button = widgets.Button(description="Интерполировать")
        self.save_button = widgets.Button(description="Сохранить координаты")
        self.add_point_button = widgets.Button(description="Добавить точки")

        self.x_input = widgets.FloatText(description="X:", step=0.1)
        self.y_input = widgets.FloatText(description="Y:", step=0.1)

        self.num_points_input = widgets.IntText(description="Кол-во точек:", value=10, min=2)
        self.step_input = widgets.FloatText(description="Шаг:", value=0.0)

        self.clear_button.on_click(self.clear_points)
        self.interpolate_button.on_click(self.interpolate_points)
        self.save_button.on_click(self.save_coordinates)

        self.add_point_button.on_click(self.add_point_manually)
        self.last_click_time = 0  
        self.double_click_delay = 0.3  
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)

    def show(self):
        display(
            widgets.VBox([
                self.fig.canvas,
                widgets.HBox([
                    self.clear_button,
                    self.interpolate_button,
                    self.save_button
                ]),
                self.coordinates_output,
                self.output
            ]),
            widgets.HBox([
                self.x_input,
                self.y_input,
                self.add_point_button,
                self.num_points_input,
                self.step_input
            ])
        )

    def save_coordinates(self, button):
        with self.coordinates_output:
            print("\nВсе координаты точек:")
            for i, (x, y) in enumerate(self.points, 1):
                print(f"Точка {i}: x={x:.2f}, y={y:.2f}")
            print("")

    def on_click(self, event):
        if event.inaxes != self.ax:
            return
        
        current_time = time.time()
        
        if (current_time - self.last_click_time) < self.double_click_delay:
            self.add_point(event.xdata, event.ydata)
            self.last_click_time = 0  
        else:
            self.last_click_time = current_time

    def add_point(self, x, y):  
        self.points.append((x, y))
        point = self.ax.plot(x, y, 'ro', markersize=5)[0]
        self.point_artists.append(point)
    
        with self.coordinates_output:
            print(f"Добавлена точка: x={x:.2f}, y={y:.2f}")
    
        self.fig.canvas.draw()

    def add_point_manually(self, button):
        x = self.x_input.value
        y = self.y_input.value
        self.add_point(x, y)
 

    def clear_points(self, button):
        for artist in self.point_artists:
            artist.remove()
        
        if self.line_artist is not None:
            self.line_artist.remove()
        
        self.points.clear()
        self.point_artists.clear()
        self.line_artist = None
        
        self.coordinates_output.clear_output()
        self.output.clear_output()

        current_clim = self.img.get_clim()
        self.img.set_array(self.image_data)
        self.img.set_clim(current_clim)
        
        self.fig.canvas.draw_idle()

    def interpolate_points(self, button):
        if len(self.points) < 3:
            with self.output:
                print("Нужно хотя бы 3 точки для интерполяции!")
            return

        with self.output:
            print(f"\nИнтерполяция выполнена через точки:")
            for i, (x, y) in enumerate(self.points, 1):
                print(f"Точка {i}: ({x:.2f}, {y:.2f})")

        points = np.array(self.points)
        x, y = points[:, 0], points[:, 1]

        t = np.zeros(len(x))
        for i in range(1, len(x)):
            dx = x[i] - x[i - 1]
            dy = y[i] - y[i - 1]
            t[i] = t[i - 1] + np.sqrt(dx ** 2 + dy ** 2)

        fx = interp1d(t, x, kind='cubic')
        fy = interp1d(t, y, kind='cubic')

        step = self.step_input.value
        n_points = self.num_points_input.value

        if step > 0:
            t_new = np.arange(t[0], t[-1], step)
        else:
            t_new = np.linspace(t[0], t[-1], n_points)

        x_new = fx(t_new)
        y_new = fy(t_new)

        if self.line_artist is not None:
            self.line_artist.remove()

        self.line_artist, = self.ax.plot(x_new, y_new, 'b-', linewidth=2)

        for artist in self.point_artists:
            artist.remove()
        self.point_artists = []

        for x_i, y_i in zip(x_new, y_new):
            pt = self.ax.plot(x_i, y_i, 'go', markersize=4)[0]
            self.point_artists.append(pt)

        with self.coordinates_output:
            print("\nКоординаты точек вдоль интерполированной линии:")
            for i, (x_i, y_i) in enumerate(zip(x_new, y_new), 1):
                print(f"{i}: x={x_i:.2f}, y={y_i:.2f}")

        self.fig.canvas.draw()