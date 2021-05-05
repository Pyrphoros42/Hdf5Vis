import numpy as np
import matplotlib.pyplot as plt


def gauss():

    len = 10000
    x = []
    # len values normally distributed
    for i in range(len):
        x.append(np.random.randn())

    plt.figure()
    plt.subplot(211)

    # int + 1 edges drawn
    hx, hy, _ = plt.hist(x=[x,2*x], bins=int(len**0.5)*3)

    # making it more beautiful
    plt.title('My first normal distribution')
    plt.grid()

    plt.subplot(212)
    plt.plot([1,2,3,4],[2,4,6,8])
    plt.savefig('first double.png')
    # show it
    plt.show()
