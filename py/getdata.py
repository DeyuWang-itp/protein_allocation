import matplotlib.pyplot as plt
import matplotlib.image as mpimg

def getdata(imageFile):
    img = mpimg.imread(imageFile)
    fig, ax = plt.subplots()
    ax.imshow(img)
    clicked_coordinates = []
    axis_point = []
    axis_label = ['x1','x2','y1','y2']
    data_point = []
    print("please selected axis point as x1 x2 y1 y2")
    def _get_pixel(event):
        x,y = event.xdata, event.ydata
        if x is not None and y is not None:
            if len(axis_point) < 4:
                ax.plot(x,y,'ro', label = 'axis label')
                axis_point.append((x,y))
            else:
                ax.plot(x,y,'go', label = 'data label')
                data_point.append((x,y))
            plt.draw()
    cid = fig.canvas.mpl_connect('button_press_event', _get_pixel)
    plt.show()
    x1 = float(input('x1'))
    x2 = float(input('x2'))
    y1 = float(input('y1'))
    y2 = float(input('y2'))
    dx = x2 - x1
    dy = y2 - y1
    x10 = axis_point[0][0]
    x20 = axis_point[1][0]
    y10 = axis_point[2][1]
    y20 = axis_point[3][1]
    # breakpoint()
    data = [((xi - x10)/(x20 - x10) * dx + x1, ((yi - y10)/(y20 - y10) * dy + y1)) for xi, yi in data_point]
    return data



