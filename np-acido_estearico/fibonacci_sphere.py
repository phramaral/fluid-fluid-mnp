import math
import numpy as np

def fibonacci_sphere(samples,r):

    points = []
    phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = math.sqrt(1 - y * y)  # radius at y
        
        theta = phi * i  # golden angle increment

        x = math.cos(theta) * radius*r
        z = math.sin(theta) * radius*r

        points.append((x, y*r, z))

    return points

n=154
r=17.50
points = np.array(fibonacci_sphere(n, r))

a = np.array(['ap']*n).reshape(n, 1)

b = np.concatenate([a, points], axis=1)

