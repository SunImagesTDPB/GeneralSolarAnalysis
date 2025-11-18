import tkinter as tk

class SimpleConverter:
    def __init__(self, root, center_x, center_y, center_tx, center_ty):
        self.root = root
        self.root.title("Калькулятор координат (без FITS)")

        self.scale = 0.6  
        self.center_x = center_x
        self.center_y = center_y
        self.center_tx = center_tx
        self.center_ty = center_ty

        tk.Label(root, text="Солнечные координаты (arcsec):").pack()
        self.tx_entry = tk.Entry(root)
        self.tx_entry.pack()
        self.ty_entry = tk.Entry(root)
        self.ty_entry.pack()
        tk.Button(root, text="Перевод в пиксели", command=self.to_pixel).pack(pady=5)

        tk.Label(root, text="Пиксельные координаты:").pack()
        self.xpix_entry = tk.Entry(root)
        self.xpix_entry.pack()
        self.ypix_entry = tk.Entry(root)
        self.ypix_entry.pack()
        tk.Button(root, text="Перевод в солнце", command=self.to_solar).pack(pady=5)

        self.result = tk.Label(root, text="Результат:")
        self.result.pack(pady=10)

    def to_pixel(self):
        try:
            tx = float(self.tx_entry.get())
            ty = float(self.ty_entry.get())
            x = (tx - self.center_tx) / self.scale + self.center_x
            y = (ty - self.center_ty) / self.scale + self.center_y
            self.result.config(text=f"Пиксели: x = {x:.1f}, y = {y:.1f}")
        except Exception as e:
            self.result.config(text=f"Ошибка: {e}")

    def to_solar(self):
        try:
            x = float(self.xpix_entry.get())
            y = float(self.ypix_entry.get())
            tx = (x - self.center_x) * self.scale + self.center_tx
            ty = (y - self.center_y) * self.scale + self.center_ty
            self.result.config(text=f"Солнце: Tx = {tx:.1f}″, Ty = {ty:.1f}″")
        except Exception as e:
            self.result.config(text=f"Ошибка: {e}")

if __name__ == "__main__":
    root = tk.Tk()
    app = SimpleConverter(root, 835//2, 502//2, -268.8, -268.8)
    root.geometry("400x350")
    root.mainloop()
