from Score import Score
from GPR import GPR
from MPR import MPR
from TS import TS
import sys

def main():
    filename = sys.argv[1]

    score = Score(filename)

    gpr = GPR(score)
    gpr.apply_rules()
    G_nodes = gpr.get_result()

    mpr = MPR(score,G_nodes)
    mpr.apply_rules()
    M_nodes = mpr.get_result()

    ts = TS(score,G_nodes,M_nodes)
    ts.apply_rules()
    T_nodes = ts.get_result()

    print_TS(T_nodes,0)

    

def print_TS(node,n):
        if node == None:
            return
        if n < 0:
            n *= -1
            for i in range(n):
                print("--", end = "")
        else:
            for i in range(n):
                print("==", end = "")
        print(node.L_end,node.dot)
        print_TS(node.primary,abs(n+1) ) 
        print_TS(node.secondary,abs(n+1) * -1)

if __name__ == "__main__":
    main()
