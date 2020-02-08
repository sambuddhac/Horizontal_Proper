import matplotlib.pyplot as plt

with open('UpperBoundIterates.txt') as f:
    lines = f.readlines()
    x = [int(line.split()[0]) for line in lines]
    y = [int(line.split()[1]) for line in lines]

fig = plt.figure(figsize=(15,10))
ax1 = fig.add_subplot(111)
ax1.set_title("SMP graph")
ax1.set_xlabel('hour')
ax1.set_ylabel('smp')
ax1.bar(x,y, width=0.7)

fig1=plt.gcf()
plt.show()
plt.draw()
