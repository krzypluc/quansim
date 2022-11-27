import matplotlib.pyplot as plt

def ComplexPlot(x, y):
    plt.plot(x, y.real, x, y.imag)
    plt.show()

def ComplexTwoPlots(x, y1, y2, title1="", title2="", title="", fname="", xmin=None, xmax=None, ymin=None, ymax=None):
    fig, ax = plt.subplots(1, 2)
    fig.set_size_inches(18.5, 10.5)
    fig.suptitle(title)
    ax[0].plot(x, y1.real, x, y1.imag)
    ax[0].title.set_text(title1)
    ax[0].set(xlim=(xmin, xmax), ylim=(ymin, ymax))
    ax[1].plot(x, y2.real, x, y2.imag)
    ax[1].title.set_text(title2)
    ax[1].set(xlim=(xmin, xmax), ylim=(ymin, ymax))
    plt.show()

def ComplexThreePlots(x, y1, y2, y3, title1="", title2="", title3="", xmin=None, xmax=None, ymin=None, ymax=None):
    fig, ax = plt.subplots(1, 3)
    fig.set_size_inches(18.5, 10.5)
    ax[0].plot(x, y1.real, x, y1.imag)
    ax[0].title.set_text(title1)
    ax[1].plot(x, y2.real, x, y2.imag)
    ax[1].title.set_text(title2)
    ax[1].set(xlim=(xmin, xmax), ylim=(ymin, ymax))
    ax[2].plot(x, y3.real, x, y3.imag)
    ax[2].title.set_text(title3)
    plt.show()
