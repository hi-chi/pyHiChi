import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rcParams.update({"font.size" : 12})

maxIter = 4

factor = 0.5
NxFull = int(factor*320)
Ny = int(factor*256)
Nz = int(factor*256)

print(NxFull, Ny, Nz)

L = 2
lenX = 40

def fn(d): return int(L*d*NxFull/lenX) * Ny * Nz
def fnlogn(d): return fn(d)*np.log(fn(d))

n_threads_arr = [1, 4, 8, 18, 28, 36]

D_ARR = range(2, 22, 2)
TIME_RTOC_ARR = {}
TIME_CTOR_ARR = {}
TIME_SOLV_ARR = {}
TIME_INIT_ARR = {}
TIME_ITER_ARR = {}
ALL_TIME_ARR = {}

for n_threads in n_threads_arr:

    D_arr = []
    time_RtoC_av = []
    time_SOLV_av = []
    time_CtoR_av = []
    time_init = []
    time_iter_av = []
    
    all_time = 0.0
    
    with open("output_%d.txt" % (n_threads), "r") as file:
        lines = [[el for el in line.split() if not el == ""] for line in file.readlines() if not line == "\n"]
        lines = [lines[i] for i in range(len(lines)-1)]
        if lines[0][0] == "1": lines = [lines[i] for i in range(2, len(lines))]
         
        STATUS = 0
        
        RtoC_sum = 0.0
        SOLV_sum = 0.0
        CtoR_sum = 0.0
        time_iter_sum = 0.0
        
        iter = -1
        
        def process_final():
            global RtoC_sum, SOLV_sum, CtoR_sum, time_iter_sum
            time_RtoC_av.append(RtoC_sum/maxIter)
            time_SOLV_av.append(SOLV_sum/maxIter)
            time_CtoR_av.append(CtoR_sum/maxIter)
            time_iter_av.append(time_iter_sum/maxIter)
            RtoC_sum = 0.0
            SOLV_sum = 0.0
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
                    all_time += float(line[5])
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
                    
            elif STATUS == 3:  # read SOLV
                if len(line) == 3 and line[1] == "PSATD:":
                    SOLV_sum += float(line[2])
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
                    all_time += float(line[6])
                    STATUS = 0 if iter == maxIter else 2
                else:
                    print("Error, status %d" % 5)
                    print(line)
                    break
                    
        process_final()
        
    time_RtoC_av = np.array(time_RtoC_av)/1000
    time_CtoR_av = np.array(time_CtoR_av)/1000
    time_SOLV_av = np.array(time_SOLV_av)/1000
    time_iter_av = np.array(time_iter_av)/1000
    time_init = np.array(time_init)/1000
        
    TIME_RTOC_ARR[n_threads] = time_RtoC_av
    TIME_CTOR_ARR[n_threads] = time_CtoR_av
    TIME_SOLV_ARR[n_threads] = time_SOLV_av
    TIME_INIT_ARR[n_threads] = time_init
    TIME_ITER_ARR[n_threads] = time_iter_av
    ALL_TIME_ARR[n_threads] = all_time
    
    print(n_threads, "threads: ALL TIME =", all_time/1000/60/60, "h =", all_time/1000/60, "min =", all_time/1000, "s")
        

    # plot time
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    ax.plot(D_arr, time_RtoC_av,  "-->", label = "RtoC")
    ax.plot(D_arr, time_CtoR_av,  "-->", label = "CtoR")
    ax.plot(D_arr, time_SOLV_av, "-->", label = "PSATD")
    ax.plot(D_arr, time_iter_av,  "-->", label = "time of one iteration")
    ax.plot(D_arr, time_init,     "-->", label = "initialization time")
    
    ax.set_xlabel("$D/L$")
    ax.set_ylabel("time [s]")
    ax.set_title("Time, number of threads = %d" % n_threads)
    ax.grid()
    ax.legend()
    
    plt.tight_layout()
    plt.savefig("time_threads_%d" % (n_threads), dpi = 500)
    plt.close(fig=fig)
    
    # plot const
    
    fnv = np.vectorize(fn)
    fnlognv = np.vectorize(fnlogn)
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    ax.plot(D_arr, time_RtoC_av/fnlognv(D_arr), "-->", label = "RtoC")
    ax.plot(D_arr, time_CtoR_av/fnlognv(D_arr), "-->", label = "CtoR")
    ax.plot(D_arr, time_SOLV_av/fnv(D_arr), "-->", label = "PSATD")
    
    ax.set_xlabel("$D/L$")
    ax.set_ylabel("time [s] divided by complexity")
    ax.set_title("Time divided by complexity (const)\n number of threads = %d" % n_threads)
    ax.grid()
    ax.legend()
    
    plt.tight_layout()
    plt.savefig("const_threads_%d" % (n_threads), dpi = 500)
    plt.close(fig=fig)
    
    
# compute scalability and create transposed arrays

RTOC_RESULTS = [[TIME_RTOC_ARR[nthr][i] for nthr in n_threads_arr] for i in range(len(D_ARR))]
SOLV_RESULTS = [[TIME_SOLV_ARR[nthr][i] for nthr in n_threads_arr] for i in range(len(D_ARR))]
CTOR_RESULTS = [[TIME_CTOR_ARR[nthr][i] for nthr in n_threads_arr] for i in range(len(D_ARR))]
INIT_RESULTS = [[TIME_INIT_ARR[nthr][i] for nthr in n_threads_arr] for i in range(len(D_ARR))]
ITER_RESULTS = [[TIME_ITER_ARR[nthr][i] for nthr in n_threads_arr] for i in range(len(D_ARR))]

def compEff(time, time1, nthr):
    return time1/(time * nthr)

SCAL_RTOC_RESULTS = [[compEff(TIME_RTOC_ARR[nthr][i], TIME_RTOC_ARR[1][i], nthr) for nthr in n_threads_arr] for i in range(len(D_ARR))]
SCAL_SOLV_RESULTS = [[compEff(TIME_SOLV_ARR[nthr][i], TIME_SOLV_ARR[1][i], nthr) for nthr in n_threads_arr] for i in range(len(D_ARR))]
SCAL_CTOR_RESULTS = [[compEff(TIME_CTOR_ARR[nthr][i], TIME_CTOR_ARR[1][i], nthr) for nthr in n_threads_arr] for i in range(len(D_ARR))]
SCAL_INIT_RESULTS = [[compEff(TIME_INIT_ARR[nthr][i], TIME_INIT_ARR[1][i], nthr) for nthr in n_threads_arr] for i in range(len(D_ARR))]
SCAL_ITER_RESULTS = [[compEff(TIME_ITER_ARR[nthr][i], TIME_ITER_ARR[1][i], nthr) for nthr in n_threads_arr] for i in range(len(D_ARR))]

# plot scalability

for i in range(len(D_ARR)):

    # plot time 
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    ax.plot(n_threads_arr, RTOC_RESULTS[i], "-->", label = "RtoC")
    ax.plot(n_threads_arr, CTOR_RESULTS[i], "-->", label = "CtoR")
    ax.plot(n_threads_arr, SOLV_RESULTS[i], "-->", label = "PSATD")
    ax.plot(n_threads_arr, ITER_RESULTS[i], "-->", label = "time of one iteration")
    ax.plot(n_threads_arr, INIT_RESULTS[i], "-->", label = "initialization time")
    
    ax.set_xlabel("number of threads")
    ax.set_ylabel("time [s]")
    ax.set_title("Time, D/L = %d" % D_ARR[i])
    ax.grid()
    ax.legend()
    
    plt.tight_layout()
    plt.savefig("time_DL_%d" % (D_ARR[i]), dpi = 500)
    plt.close(fig=fig)
    
    #plot scalability
       
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    ax.plot(n_threads_arr, SCAL_RTOC_RESULTS[i], "-->", label = "RtoC")
    ax.plot(n_threads_arr, SCAL_CTOR_RESULTS[i], "-->", label = "CtoR")
    ax.plot(n_threads_arr, SCAL_SOLV_RESULTS[i], "-->", label = "PSATD")
    ax.plot(n_threads_arr, SCAL_ITER_RESULTS[i], "-->", label = "time of one iteration")
    ax.plot(n_threads_arr, SCAL_INIT_RESULTS[i], "-->", label = "initialization time")
    
    ax.set_xlabel("number of threads")
    ax.set_ylabel("efficiency")
    ax.set_title("Scalability, D/L = %d" % D_ARR[i])
    ax.set_ylim((0.0, 1.2))
    ax.grid()
    ax.legend()
    
    plt.tight_layout()
    plt.savefig("scal_DL_%d" % (D_ARR[i]), dpi = 500)
    plt.close(fig=fig)


# create csv file

def writeARRToFile(nameFile, arr_of_res_arr, label_table):
    with open(nameFile, "w") as file:
        for res_arr, label in zip(arr_of_res_arr,
                                  ["RtoC", "PSATD", "CtoR", "initialization", "one iteration"]
                                 ):
            file.write("%s [s] of %s\n" % (label_table, label))
            file.write(";;number of threads\n")
            file.write("D/L;grid size;")
            for nthr in n_threads_arr:
                file.write("%d;" % nthr)
            file.write("\n")
            for i in range(len(D_ARR)):
                d = D_ARR[i]
                file.write("%d;%dx%dx%d;" % (d, int(L*d*NxFull/lenX), Ny, Nz))
                for j in range(len(n_threads_arr)):
                    file.write(("%0.2f;" % res_arr[i][j]).replace(".", ","))
                file.write("\n")
            file.write("\n\n")
    
writeARRToFile("time_results.csv", [RTOC_RESULTS, SOLV_RESULTS, CTOR_RESULTS, INIT_RESULTS, ITER_RESULTS], "TIME")   
writeARRToFile("scal_results.csv", [SCAL_RTOC_RESULTS, SCAL_SOLV_RESULTS, SCAL_CTOR_RESULTS, SCAL_INIT_RESULTS, SCAL_ITER_RESULTS], "SCALABILITY")   