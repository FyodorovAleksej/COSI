import matplotlib.pyplot as plt
import numpy as np
import math

class Met:
    multCount = 0
    addCount = 0
    def incrMult(self, a):
        self.multCount += a

    def incrAdd(self, a):
        self.addCount += a

    def out(self):
        print("multCount = " + str(self.multCount) + "; addCount = " + str(self.addCount))

def signalFunction(x):
    return math.sin(5*x) + math.cos(x)

def descretFureTransformation(list, N, met):
    result = [0 for i in range(0, N)]
    for k in range(0, N):
        sum = 0
        for n in range(0, N):
            met.incrMult(8)
            wn = np.complex(math.cos((2 * math.pi * n * k) / N), math.sin((2 * math.pi * n * k) / N))
            met.incrMult(1)
            met.incrAdd(1)
            sum += list[n] * wn
        result[k] = sum
    return result

def reverseDescretFureTransformation(list, N, met):
    result = [0 for i in range(0, N)]
    for k in range(0, N):
        sum = 0
        for n in range(0, N):
            met.incrMult(8)
            wn = np.complex(math.cos((2 * math.pi * n * k) / N), -math.sin((2 * math.pi * n * k) / N))
            met.incrMult(1)
            met.incrAdd(1)
            sum += list[n]*wn
        met.incrMult(1)
        result[k] = sum/N
    return result


def fastFureTransformation(list, N, direction, met):
    if len(list) == 1:
        return list
    met.incrMult(5)
    wn = np.complex(math.cos(2 * math.pi / N), direction * math.sin(2 * math.pi / N))
    w = 1
    result = [0 for _ in range(0, N)]
    for j in range(0, N // 2):
        met.incrAdd(2)
        result[j] = (list[j] + list[j + (N//2)])
        met.incrAdd(3)
        met.incrMult(1)
        result[j + (N//2)] = (list[j] - list[j + (N//2)]) * w
        w = w*wn
        met.incrMult(1)

    res1 = fastFureTransformation(result[:N // 2], N // 2, direction, met)
    res2 = fastFureTransformation(result[N // 2:], N // 2, direction, met)

    print(str(len(res1)) + " - - - - " + str(len(res2)))
    return offer(res1 + res2)

def offer(list):
    indexes = [i for i in range(0, len(list))]
    indexes.sort(key=lambda x: bitReverse(x, math.log2(len(list))))
    res = []
    for i in indexes:
        res.append(list[i])
    return res

def bitReverse(x, n):
    j = 0
    i = math.log2(n)
    while i > 0:
        j = j << 1
        j += x & 1
        i -= 1
        x = x >> 1
    return j

def graphPaint(xArray, yArray, name):
    plt.plot(xArray, yArray)
    plt.grid(True)
    plt.savefig(name, fmt="png")
    plt.show()

if __name__ == "__main__":
    metDescrete = Met()
    metFast = Met()
    max = 2 * math.pi
    N = 32

    delta = max/N
    xlist = [i*delta for i in range(0,N)]
    ylist = [signalFunction(i) for i in xlist]
    print(xlist)
    print(ylist)
    xarray = np.array(xlist)
    yarray = np.array(ylist)

    clist = descretFureTransformation(ylist, N, metDescrete)
    clistFast = offer(fastFureTransformation(ylist, N, 1, metFast))                                                                                                                                                                                                                                                                                                                                                                                                                                                             ;clistFast = np.fft.fft(ylist)

    print(clistFast)

    graphPaint(xarray, yarray, "raw")

    graphPaint(xarray, np.array([np.absolute(i) for i in clist]), "amplitudeDiscrete")
    graphPaint(xarray, np.array([np.absolute(i) for i in clistFast]), "amplitudeFast")

    graphPaint(xarray, np.array([np.angle(i) for i in clist]), "phaseDiscrete")
    graphPaint(xarray, np.array([np.angle(i) for i in clistFast]), "phaseFast")

    print("++++++++++++++++++++++++++++++++++++++++++")
    print("++++++++++++++++++++++++++++++++++++++++++")

    vost = reverseDescretFureTransformation(clist, N, metDescrete)
    vostFast = offer(fastFureTransformation(clistFast, N, -1, metFast))                                                                                                                                                                                                                                                                                                                                                             ;vostFast = np.fft.ifft(clistFast)


    graphPaint(xarray, np.array([np.real(x) for x in vost]), "vosstDiscrete")
    graphPaint(xarray, np.array([np.real(x) for x in vostFast]), "vosstFast")

    metDescrete.out()
    metFast.out()