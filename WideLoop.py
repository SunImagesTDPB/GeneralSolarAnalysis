from InteractiveLoopTracer import *

class WideLoop(InteractiveLoopTracer):
    def __init__(self, image_data):
        super().__init__(image_data)
        
        self.line_artists = []  
        
        self.clear_all_button = widgets.Button(description="Очистить изображение")
        self.clear_all_button.on_click(self.clear_all)

        self.create_line_button = widgets.Button(description="Добавить полосы")
        self.create_line_button.on_click(self.create_line)

        self.lenght_line_input = widgets.FloatText(description="Длина полосы:", value=5.0, layout=widgets.Layout(width='150px'))

    def show(self):
        display(
            widgets.VBox([
                self.fig.canvas,
                widgets.HBox([
                    self.clear_all_button,
                    self.interpolate_button,
                    self.create_line_button
                ]),
                self.coordinates_output,
                self.output
            ]),
            widgets.HBox([
                self.x_input,
                self.y_input,
                self.add_point_button,
                self.num_points_input,
                self.step_input,
                self.lenght_line_input
            ])
        )

    def clear_all(self, button):
        super().clear_points(button)
        
        for artist in self.line_artists:
            artist.remove()
        self.line_artists = []
        
        self.fig.canvas.draw_idle()

    def create_line(self, button):

        length = self.lenght_line_input.value

        x_line = self.line_artist.get_xdata()
        y_line = self.line_artist.get_ydata()

        dx = np.gradient(x_line)
        dy = np.gradient(y_line)
        
        norm = np.sqrt(dx**2 + dy**2)
        nx = -dy / norm  
        ny = dx / norm

        for sign in [-1, 1]:
            x_ends = x_line + sign * length * nx
            y_ends = y_line + sign * length * ny

            for i in range(len(x_line)):
                line = self.ax.plot(
                    [x_line[i], x_ends[i]], 
                    [y_line[i], y_ends[i]], 
                    color='red', 
                    linewidth=1.2,
                    alpha=0.6
                )[0]
                self.line_artists.append(line)

        self.fig.canvas.draw()