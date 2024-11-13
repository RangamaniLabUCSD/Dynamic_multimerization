import Main_Analyze_Trajectory as main
import sys
N = int(sys.argv[1])
Rval = int(sys.argv[2])
repid = int(sys.argv[3])
print(len(sys.argv), flush = True)
main.frontend_cross_DynamicMultimer_set_Rval_repid(N, Rval, repid, sys.argv[4], sys.argv[5])