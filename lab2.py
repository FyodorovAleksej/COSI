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


def correlation(yarray1, yarray2, N, met):
    result = [0 for i in range(0 , N)]
    for m in range(0 ,N):
        for h in range(0, N):
            if (m + h < N):
                result[m] += yarray1[h] * yarray2[m + h] #CORRELATION
                met.incrAdd(2)
                met.incrMult(1)
            else:
                result[m] += yarray1[h] * yarray2[m + h - N]
                met.incrAdd(2)
                met.incrMult(1)
        result[m] /= N
        met.incrMult(1)
    return result


def convolution(yarray1, yarray2, N, met):
    result = [0 for i in range(0, N)]
    for m in range(0, N):
        for h in range(0, N):
            if (m - h >= 0):
                result[m] += yarray1[h] * yarray2[m - h]
                met.incrAdd(2)
                met.incrMult(1)
            else:
                result[m] += yarray1[h] * yarray2[m - h + N]
                met.incrAdd(2)
                met.incrMult(1)
        result[m] /= N
        met.incrMult(1)
    return result

def yFunction(x):
    return math.sin(5*x)

def zFunction(x):
    return math.cos(x)


def fastFureTransformation(list, N, direction, met):
    if len(list) == 1:
        return list
    wn = np.complex(math.cos(2 * math.pi / N), direction * math.sin(2 * math.pi / N))
    w = 1
    result = [0 for _ in range(0, N)]
    for j in range(0, N // 2):
        result[j] = (list[j] + list[j + (N//2)])
        met.incrAdd(2)
        met.incrMult(1)
        result[j + (N//2)] = (list[j] - list[j + (N//2)]) * w
        met.incrAdd(2)
        met.incrMult(1)
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
    metFast = Met()
    metCorrel = Met()
    metConvol = Met()
    max = 2 * math.pi
    N = 16

    delta = max / N
    xlist1 = [i * delta for i in range(0, N)]
    ylist1 = [yFunction(i) for i in xlist1]

    xlist2 = [i * delta for i in range(0, N)]
    ylist2 = [zFunction(i) for i in xlist2]

    print(xlist1)
    print(ylist1)
    print("\n\n")
    print(xlist2)
    print(ylist2)

    xarray1 = np.array(xlist1)
    yarray1 = np.array(ylist1)

    xarray2 = np.array(xlist2)
    yarray2 = np.array(ylist2)

    clistFast1 = offer(fastFureTransformation(ylist1, N, 1, metFast))
    clistFast1 = np.fft.fft(ylist1)
    clistFast2 = offer(fastFureTransformation(ylist2, N, 1, metFast))
    clistFast2 = np.fft.fft(ylist2)

    correl = np.array(correlation(yarray1, yarray2, N, metCorrel))
    convol = np.array(convolution(yarray1, yarray2, N, metConvol))

    plt.plot(xarray1, yarray1)
    plt.plot(xarray2, yarray2)
    plt.grid(True)
    plt.savefig("raw", fmt="png")
    plt.show()

    graphPaint(xarray1, convol, "convolution")
    graphPaint(xarray1, correl, "correlation")

    resultBPF = [clistFast1[i] * clistFast2[i] for i in range(0, N)]
    vostFast = offer(fastFureTransformation(resultBPF, N, -1, metFast));
    vostFast = np.fft.ifft(resultBPF)
    graphPaint(xarray1, np.array([np.real(x) for x in vostFast]), "vosstFast")
    print("Convolution operations:")
    metConvol.out()
    print("Correlation operations:")
    metCorrel.out()
    print("FFT opers:")
    metFast.out()
