import matplotlib
import matplotlib.pyplot as plt

maxIter = 32

D_arr = []
time_RtoC_av = []
time_PSATD_av = []
time_CtoR_av = []
time_init = []
time_iter_av = []

with open("restime.txt", "r") as file:
    lines = [[el for el in line.split() if not el == ""] for line in file.readlines() if not line == "\n"]
     
    STATUS = 0
    
    RtoC_sum = 0.0
    PSATD_sum = 0.0
    CtoR_sum = 0.0
    time_iter_sum = 0.0
    
    iter = -1
    
    def process_final():
        global RtoC_sum, PSATD_sum, CtoR_sum, time_iter_sum
        time_RtoC_av.append(RtoC_sum/maxIter)
        time_PSATD_av.append(PSATD_sum/maxIter)
        time_CtoR_av.append(CtoR_sum/maxIter)
        time_iter_av.append(time_iter_sum/maxIter)
        RtoC_sum = 0.0
        PSATD_sum = 0.0
        CtoR_sum = 0.0
        time_iter_sum = 0.0
        
    
    for line in lines:
        
        if STATUS == 0:  # read D
            if len(line) == 4 and line[0] == "D":
                if (not iter == -1):
                    process_final()
                iter = 0
                D_arr.append(float(line[2]))
                STATUS = 1
            else:
                print("Error, status %d" % 0)
                print(line)
                break
        
        elif STATUS == 1:  # read init
            if len(line) == 7 and line[3] == "init":
                time_init.append(float(line[5]))
                STATUS = 2
            else:
                print("Error, status %d" % 1)
                print(line)
                break
                
        elif STATUS == 2:  # read RtoC
            if len(line) == 3 and line[1] == "RtoC:":
                RtoC_sum += float(line[2])
                STATUS = 3
            else:
                print("Error, status %d" % 2)
                print(line)
                break
                
        elif STATUS == 3:  # read PSATD
            if len(line) == 3 and line[1] == "PSATD:":
                PSATD_sum += float(line[2])
                STATUS = 4
            else:
                print("Error, status %d" % 3)
                print(line)
                break

        elif STATUS == 4:  # read CtoR
            if len(line) == 3 and line[1] == "CtoR:":
                CtoR_sum += float(line[2])
                STATUS = 5
            else:
                print("Error, status %d" % 4)
                print(line)
                break    

        elif STATUS == 5:  # read time iter
            if len(line) == 8 and line[4] == "iter" and int(line[3]) == iter:
                iter += 1
                time_iter_sum += float(line[6])
                STATUS = 0 if iter == maxIter else 2
            else:
                print("Error, status %d" % 5)
                print(line)
                break
                
    process_final()
    
    
print("D            ", D_arr        )
print("time_RtoC_av ", time_RtoC_av )
print("time_PSATD_av", time_PSATD_av)
print("time_CtoR_av ", time_CtoR_av )
print("time_init    ", time_init    )
print("time_iter_av ", time_iter_av )


fig = plt.figure()
ax = fig.add_subplot(1,1,1)

ax.plot(D_arr, time_RtoC_av, "-->", label = "RtoC")
ax.plot(D_arr, time_CtoR_av, "-->", label = "CtoR")
ax.plot(D_arr, time_PSATD_av, "-->", label = "PSATD")
ax.plot(D_arr, time_iter_av, "-->", label = "time of one iteration")
ax.plot(D_arr, time_init, "-->", label = "initialization time")

ax.set_xlabel("$D/L$")
ax.set_ylabel("time, ms")

ax.legend()
plt.savefig("graph", dpi = 500)

    