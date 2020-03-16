import statistics
import matplotlib.pyplot as plt

file = open("GA_30_Sphere.txt","r")
dateGA = list()
for i in range(0,30):
    number = file.readline().split("\n")[0]
    dateGA.append(float(number))

print(dateGA)

file1 = open("HCB_30_Sphere.txt","r")
dateHCB = list()
for i in range(0,30):
    number = file1.readline().split("\n")[0]
    dateHCB.append(float(number))

print(dateHCB)

file2 = open("HCF_30_Sphere.txt","r")
dateHCF = list()
for i in range(0,30):
    number = file2.readline().split("\n")[0]
    dateHCF.append(float(number))

print(dateHCF)

file3 = open("SA_30_Sphere.txt","r")
dateSA = list()
for i in range(0,30):
    number = file3.readline().split("\n")[0]
    dateSA.append(float(number))

print(dateSA)
run = list(range(1, 31))
print(run)

plt.plot(run, dateGA, color='r')
plt.plot(run, dateHCB, color='orange')
plt.plot(run, dateHCF, color='b')
plt.plot(run, dateSA, color='g')
plt.xlabel('Running Times')
plt.ylabel('Minum Returned')
plt.title('Evolution of minimum on Sphere 30 D')
plt.show()

#print(str(round(min(date),5)) + " & " + str(round(max(date),5)) + " & " + str(round(statistics.mean(date), 5)) + " & " + str(round(statistics.median(date), 5)) + " & " + str(round(statistics.pstdev(date), 5)))
