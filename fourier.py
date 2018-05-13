# Load and analyse via a fourier transform on data that looks like:
# === T(ime) Y   Z
# 1.0        2.0 3.0
# 0.5        1.0 1.5
# 2.1        444 5.5
# Then find peaks and finally present the transforms on a nice graph.


import numpy as np
import matplotlib.pyplot as plt
import sys

from scipy import signal
PRINT_ALL = True

# Little macro things
def pair(i):
    return [x[i], y[i]]
def paird(x, y, i):
    return [x[i], y[i]]


argc = len(sys.argv)
# get data file
if (argc < 2):
    print("Usage: %s <file-path> [min frequency threshold]" % (sys.argv[0]))
    exit(1)

PEAK_MIN = 50

# Load data
path = sys.argv[1]
if (argc > 2):
    PEAK_MIN = float(sys.argv[2])

def loadArrays(path):
    f = open(path, "r")
    FILE = f.read()
    f.close()
    # Load them into slow arrays
    tempArrays = {}
    lines = FILE.split("\n")
    cats = None
    for line in lines:
        if len(line) < 3: continue

        if line[:3] == "===":
            # Table header
            cats = (line[3:]).strip().split()
            for cat in cats:
                tempArrays[cat] = []
            continue


        if cats:
            vals = line.split()
            if len(vals) != len(cats):
                print("Error, value count mismatch on line:", line)
                exit(1)

            for i, cat in enumerate(cats):
                tempArrays[cat].append(float(vals[i]))


    # Make them numpy arrays
    for arr in tempArrays:
        tempArrays[arr] = np.array(tempArrays[arr])


    # Always assume that the first of cats is time
    return tempArrays, cats

tempArrays, cats = loadArrays(path)

def fourier(time, y):
    tStep = np.polyfit(np.arange(time.size), time, 1)[0]
    freqAxis = np.fft.rfftfreq(time.size, d=tStep)
    fft = np.abs(np.fft.rfft(y))

    return np.array([freqAxis, fft])

tStep = np.polyfit(np.arange(tempArrays[cats[0]].size), tempArrays[cats[0]], 1)[0].round(6)

tArray = tempArrays[cats[0]]
tArray = np.arange(tArray.size)*tStep + tArray[0]
tempArrays[cats[0]] = tArray


# Calculate FFT of all other arrays
FFTs = {}

for cat in cats[1:]:
    FFTs[cat] = fourier(tempArrays[cats[0]], tempArrays[cat])



# Find maxima
Peaks = {}
maxYindex = int(4 / tStep)# Assume frequency is < 1 Hz
for FFT in FFTs:
    if PRINT_ALL: print("Peak analysis for %s:" % (FFT))
    x = FFTs[FFT][0]
    y = FFTs[FFT][1]
    peaks = signal.find_peaks_cwt(y[:maxYindex], np.array([0.01]) / tStep)
    peakRemove = [] # Peaks to remove after for loop
    for i, p in enumerate(peaks):
        tMax = (y[p - 1] <= y[p] and (y[p] - y[p + 1]))
        # If sloping down right, return False
        # If sloping down left, return -ve
        # If min, return False
        # If max, return True
        # if tMax: continue
        counter = 0

        while (not tMax):
            if (y[p] == y[p + 1]): break # Otherwise flat tops will create infinite loops
            if tMax < 0: p += 1
            if tMax == False: p -= 1
            tMax = (y[p - 1] <= y[p] and (y[p] - y[p + 1]))
            counter += 1
            # if (counter > 10):
            if (counter > 10):
                if PRINT_ALL: print("-? Removing max at %d (%f)" % (peaks[i], y[peaks[i]]))
                peakRemove.append(i)
                continue

        if (y[p] < PEAK_MIN):
            if PRINT_ALL: print("-- Removing too-small at %d (%f)" % (peaks[i], y[peaks[i]]))
            peakRemove.append(i)
            continue

        if counter != 0:
            if (p < peaks[i]):
                if PRINT_ALL: print("<< ", end="")
            else:
                if PRINT_ALL: print(">> ", end="")
            if PRINT_ALL: print("Moving peak at %d (%f) to %d (%f)" % (peaks[i], y[peaks[i]], p, y[p]))
            peaks[i] = p
        else:
            if PRINT_ALL: print("== Leaving peak at %d (%f)." % (p, y[p]))

    peaks = np.delete(peaks, peakRemove)

    Peaks[FFT] = [[j.round(4) for j in pair(i)] for i in peaks]

print("Final peaks:")
for set in Peaks:
    print("   %s:" % set)
    print("   ", Peaks[set], end = "")
    print()


# Plot all transforms
i = 1
count = len(FFTs)

for FFT in FFTs:
    ax = plt.subplot(count, 1, 1)
    if (i == 1):
        ax.set_xlim(left=0, right=100)
        ax.grid()
        plt.ylabel(", ".join(cats[1:]))
    # plt.subplot(count, 2, i * 2 - 1)
    if i == 1: plt.title("Displacement vs time of: %s. File: %s" % (", ".join(cats[1:]), path))
    if i == count: plt.xlabel("Time (s)")
    y = tempArrays[FFT]
    x = np.arange(y.size) * tStep
    plt.plot(x, y)
    # plt.subplot(count, 2, i * 2)
    ax = plt.subplot(2, count, count + i)
    # if i == 1: plt.title("Fourier transforms of: %s" % (", ".join(cats[1:])))
    plt.title("Fourier transform of %s" % (FFT))
    x = FFTs[FFT][0]
    y = FFTs[FFT][1]
    plt.xlabel("Frequency (Hz)")
    ax.grid()
    ax.set_xlim(left=0, right=3)

    plt.plot(x, y)
    i += 1

# input("Enter to show graphs")
plt.show()




















##
