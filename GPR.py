import math
from os import write
import numpy as np
import xml.etree.ElementTree as et
import xml.dom.minidom as md
from GTTMRuleSet import GTTMRuleSet

class GPR_node:
    def __init__(self, L_end, rest, ioi, regi, dyn, arti, leng, id):
        self.L_end = L_end
        self.rest = rest #消音時刻から発音時刻までの間隔
        self.ioi = ioi #発音時刻間隔
        self.regi = regi #音高の差
        self.dyn = dyn #ダイナミクスの差
        self.arti = arti #楽譜上の音符の長さと実際に演奏された音の長さの比
        self.leng = leng #音価の差
        self.rule = GPR_rule() #適応ルール
        self.rule.id = id

class GPR_rule:
    def __init__(self):
        self.boundary = 0
        self.low_boundary = 0
        self.high_boundary = 0
        self.GPR = {"1":0, "2a":0, "2b":0, "3a":0, "3b":0, "3c":0, "3d":0, "4":0, "5":0, "6":0}
        self.id = ""

class Rules:
    def __init__(self,rest,ioi,regi,dyn,arti,leng):
        self.rest = rest
        self.ioi = ioi
        self.regi = regi
        self.dyn = dyn
        self.arti = arti
        self.leng = leng

    def set_nodes(self,nodes,B_low):
        self.nodes = nodes
        self.B_low = B_low
    
    #単音をグループとしない
    def GPR1(self,i):
        if self.B_low[i-1] <= self.B_low[i] and self.B_low[i+1] <= self.B_low[i] and self.nodes[i-1].rule.GPR["1"] == 0 and self.nodes[i+1].rule.GPR["1"] == 0:
            return 1
        else:
            return 0

    # GPR2:連続した４音において，各音のon-setとoff-setに注目する．
    #(a)第二音の終了時点から第三音の開始時点のまでの間
    #(b)第二音と第三音の開始時点までの間
    #が他の箇所（１と２，３と４)よりも大きければそこは境界になりやすい．
    def GPR2a(self, i):
        if i == len(self.rest)-2 or i == 1:
            return 0
        start = i-1
        end = i+1
        
        if self.rest[start] < self.rest[i] and self.rest[i] > self.rest[end]:
            return 1
        else:
            return 0
    
    def GPR2b(self, i):
        if i == len(self.rest)-2 or i == 1:
            return 0
        start = i-1
        end = i+1
        
        if self.ioi[start] < self.ioi[i] and self.ioi[i] > self.ioi[end]:
            return 1
        else:
            return 0
    
    # GPR3:連続した４音において，２と３の次の属性：
    #(a)音程(register)
    #(b)強弱(dynamics)の差
    #(c)アーティキュレーション(スタッカート，レガート)の差
    #(d)音価(length)
    #が他の箇所（１と２，３と４)よりも大きければそこは境界になりやすい．
    def GPR3a(self, i):
        if i == len(self.rest)-2 or i == 1:
            return 0
        start = i-1
        end = i+1
        
        if self.regi[start] < self.regi[i] and self.regi[i] > self.regi[end]:
            return 1
        else:
            return 0
    
    def GPR3b(self, i):
        if i == len(self.rest)-2 or i == 1:
            return 0
        start = i-1
        end = i+1

        if self.dyn[start] == 0 and self.dyn[i] != 0 and self.dyn[end] == 0:
            return 1
        else:
            return 0

    def GPR3c(self, i):
        if i == len(self.rest)-2 or i == 1:
            return 0
        start = i-1
        end = i+1

        if self.arti[start] == 0 and self.arti[i] != 0 and self.arti[end] == 0:
            return 1
        else:
            return 0

    def GPR3d(self, i):
        if i == len(self.rest)-2 or i == 1:
            return 0
        start = i-1
        end = i+1

        if self.leng[start] == 0 and self.leng[i]!= 0 and self.leng[end] == 0:
            return 1
        else:
            return 0

    #(GPR2)と(GPR3)の効果がより強く働くところが境界になりやすい．
    def GPR4(self, i, start, end):
        if i == len(self.rest)-2 or i == 1:
            return 0
        start = i-1
        end = i+1
        
        #print(self.rest[i])
        #print(i,start,end)
        if self.rest[i] != 0:
            P_rest = self.rest[i]/sum(self.rest[start:end+1])
        else :
            P_rest = 0
        if self.ioi[i] != 0:
            P_ioi = self.ioi[i]/sum(self.ioi[start:end+1])
        else:
            P_ioi = 0
        if self.dyn[i] != 0:
            P_dyn = self.dyn[i]/sum(self.dyn[start:end+1])
        else:
            P_dyn = 0
        tmp = sum(self.regi[start:end+1])
        if tmp > 0:
            P_regi = self.regi[i] / tmp
        else :
            P_regi = 0
        if self.arti[i] != 0:
            P_arti = self.arti[i]/sum(self.arti[start:end+1])
        else:
            P_arti = 0

        if max(P_rest, P_ioi, P_dyn, P_regi, P_arti) > self.T4:
            return 1
        else:
            return 0
        

    #曲は同じ長さに２分割にグループ分けされるのが望ましい．
    def GPR5(self, i, start, end):
        T_mid = sum(self.ioi[start:end])
        T_i = sum(self.ioi[start:i])
        sigma = self.adjust_sigma * T_mid

        D = (1 / math.sqrt(2 * math.pi * (sigma**2))) * math.exp((-1 * (T_i-(T_mid / 2)) ** 2 )/(2 * (sigma**2) ))

        return D
    
    #繰り返しで現れる楽句は同じグループ構造になるのが望ましい.
    def GPR6(self, i):
        D = 0
        for j in range(len(self.ioi)):
            for r in range(int(len(self.ioi)/2)):
                mij = self.GPR6_m(i,j)
                if mij == 0:
                    D += self.GPR6_G(i,j,r,True) * (1-self.Ws)
                elif mij == 1:
                    D += self.GPR6_G(i,j,r,False) * self.Ws
                elif mij == 2:
                    D += (self.GPR6_G(i,j,r,True) * (1-self.Ws)) + (self.GPR6_G(i,j,r,False) * self.Ws)
                
        return D

    def GPR6_m(self, i, j):
        qi = self.GPR6_q(i)
        qj = self.GPR6_q(j)
        if not ((i == 0 or i == len(self.ioi)-1 ) or (j == 0 or j == len(self.ioi)-1)):
            if qi != qi-self.ioi[i] and qj != qj-self.ioi[j] and qi != qi+self.ioi[i+1] and qj != qj+self.ioi[i+1]:
                mij = 0
            elif qi == qi-self.ioi[i] and qj == qj-self.ioi[j] and qi != qi+self.ioi[i+1] and qj != qj+self.ioi[j+1]:
                mij = 1
            elif qi != qi-self.ioi[i] and qj != qj-self.ioi[j] and qi == qi+self.ioi[i+1] and qj == qj+self.ioi[j+1]:
                mij = 2
            else:
                mij = 3
        else :
            mij = 3

        return mij

    def GPR6_q(self, x):
        q = 0.0
        for ioi_i in self.ioi[:x]:
            q += ioi_i
        q /= self.div
        return math.floor(q)

    def GPR6_x(self, i, r):
        x = 0
        qi = self.GPR6_q(i)
        qi_r = self.GPR6_q(i+r)
        for x in range(len(self.ioi)):
            qj = self.GPR6_q(x)
            if qi <= qj and qj <= qi_r:
                x += 1

        return x

    def GPR6_y(self, i, j, r):
        k_sum = 0
        l_sum = 0
        y = 0
        qi = self.GPR6_q(i)
        qj = self.GPR6_q(j)
        for k in self.ioi:
            k_sum += k
            for l in self.ioi:
                l_sum += l
                if abs(qi - qj) * self.div == abs(k_sum - l_sum):
                    y += 1
            l_sum = 0
        
        return y

    def GPR6_z(self, i, j, r):
        k_sum = 0
        l_sum = 0
        z = 0
        qi = self.GPR6_q(i)
        qj = self.GPR6_q(j)
        for k in self.ioi:
            k_sum += k
            for l in self.ioi:
                l_sum += l
                if abs(qi - qj) * self.div == abs(k_sum - l_sum):
                    if self.regi[i] == self.regi[j]:
                        z += 1
            l_sum = 0

        return z

    def GPR6_yz(self, i, j, r):
        k_sum = 0
        l_sum = 0
        y = 0
        z = 0
        qi = self.GPR6_q(i)
        qj = self.GPR6_q(j)
        for k in self.ioi:
            k_sum += k
            for l in self.ioi:
                l_sum += l
                if abs(qi - qj) * self.div == abs(k_sum - l_sum):
                    y += 1
                    if self.regi[i] == self.regi[j]:
                        z += 1
            l_sum = 0
        
        return y,z
                        
    def GPR6_G(self,i, j, r, is_start):
        if not is_start:
            i = i-r

        W_rdash = 1 - self.Wm
        """
        z = self.GPR6_z(i,j,r)
        y = self.GPR6_y(i,j,r)
        """
        y,z = self.GPR6_yz(i,j,r)
        x_i = self.GPR6_x(i,r)
        x_j = self.GPR6_x(j,r)
        #print(z,y,x_i,x_j)
        if y == 0:
            y = 1
            z = 0
        G = z / y * W_rdash * (r**self.Wl) + y / (x_i + x_j) * self.Wm * (r**self.Wl)

        return G

    def set_param(self,param):
        self.Wm = param["Wm"]
        self.Wl = param["Wl"]
        self.Ws = param["Ws"]
        self.div = param["div"]
        self.T4 = param["T4"]
        self.T_low = param["T_low"]
        self.adjust_sigma = param["adjust_sigma"]

class GPR(GTTMRuleSet):
    param = {"Wm":0.5,"Wl":0.5,"Ws":0.5,"S2a":0.5,"S2b":0.5,"S3a":0.5,"S3b":0.5,"S3c":0.5,"S3d":0.5,
            "S4":0.5,"S5":0.5,"S6":0.5,"T4":0.5,"T_low":0.5,"adjust_sigma":0.05,"div":1.0}

    def __init__(self,score,**kwargs):
        self.set_score(score)
        self.nodes = self.__set_GPR()
        for key,value in kwargs.items():
            self.param[key] = value

    def __set_GPR(self):
        nodes = []
        nodes.append(GPR_node(-1,0,0,0,0,0,0,self.score.get_id(0)))
        nodes[0].rule.low_boundary = 1
        nodes[0].rule.GPR["1"] = 1
        for i in range(self.score.get_score_length()-1):
            L_end = self.score.get_L_end(i)
            rest = self.score.get_rest(i)
            ioi = self.score.get_ioi(i)
            regi = self.score.get_regi(i)
            dyn = self.score.get_dyn(i) #現状指定する方法がない
            arti = self.score.get_arti(i) #同上　
            leng = self.score.get_leng(i)
            id = self.score.get_id(i+1)
            nodes.append(GPR_node(L_end, rest, ioi, regi, dyn, arti, leng, id))
        nodes.append(GPR_node(-1,0,0,0,0,0,0,""))
        nodes[-1].rule.low_boundary = 1
        nodes[-1].rule.GPR["1"] = 1

        return nodes

    def apply_rules(self):
        self.__calcGPR(self.nodes)

    def __calcGPR(self,nodes):
        rest = []
        ioi = []
        regi = []
        dyn = []
        arti = []
        leng = []
        for i in range(len(nodes)):
            rest.append(nodes[i].rest)
            ioi.append(nodes[i].ioi)
            regi.append(nodes[i].regi)
            dyn.append(nodes[i].dyn)
            arti.append(nodes[i].arti)
            leng.append(nodes[i].leng)
        
        self.rules = Rules(rest,ioi,regi,dyn,arti,leng)
        self.rules.set_param(self.param)
        rules = self.rules
        B_low = np.array([0.0] * len(nodes))
        gpr6 = np.array([0.0] * len(nodes))

        apply_rule = {"2a":rules.GPR2a, "2b":rules.GPR2b, "3a":rules.GPR3a, "3b":rules.GPR3b, "3c":rules.GPR3c, "3d":rules.GPR3d, "6":rules.GPR6}
        for i in range(1,len(nodes)-1):
            for key,rule_func in apply_rule.items():
                if key == "6":
                    gpr6[i] = rule_func(i)
                else:
                    nodes[i].rule.GPR[key] = rule_func(i)

        gpr6[1:-1] = self.min_max(gpr6[1:-1])
        for i in range(1,len(nodes)-1):
            nodes[i].rule.GPR["6"] = gpr6[i]
            for k in apply_rule.keys():
                B_low[i] += nodes[i].rule.GPR[k]*self.param["S"+k]
        B_low[0] = 1.0
        B_low[-1] = 1.0

        rules.set_nodes(nodes,B_low)
        for i in range(1,len(nodes)-1):
            nodes[i].rule.GPR["1"] = rules.GPR1(i)
            nodes[i].rule.low_boundary = 1 if nodes[i].rule.GPR["1"] == 1 and B_low[i] > self.param["T_low"] else 0
        

        self.__set_boundary(nodes, 1, 0, len(nodes))

    def __set_boundary(self,nodes,level,start,end):
        if start+1 >= end-1 :
            return 

        gpr5 = np.array([0.0] * (end-start))
        bond_number = []
        for i in range(start+1,end-1):
            if nodes[i].rule.low_boundary == 1 and nodes[i].rule.boundary == 0:
                bond_number.append(i)

        if len(bond_number) == 0:
            return
        
        B_high = np.array([0.0] * len(nodes))
        for i in range(end-start):
            gpr5[i] = self.rules.GPR5(i+start,start,end)
        gpr5 = self.min_max(gpr5)

        B_high_list = ["2a","2b","3c","3b","3c","3d","4","5","6"]
        for i in bond_number:
            nodes[i].rule.GPR["4"] = self.rules.GPR4(i,start,end)
            nodes[i].rule.GPR["5"] = gpr5[i-start]
            for k in B_high_list:
                B_high[i] += nodes[i].rule.GPR[k]*self.param["S"+k]

        B_high[bond_number] = self.min_max(B_high[bond_number])
        for i in bond_number:
            nodes[i].rule.high_boundary = nodes[i].rule.low_boundary * B_high[i]

        high_boundary = bond_number[0]
        for num in bond_number:
            if nodes[num].rule.high_boundary > nodes[high_boundary].rule.high_boundary:
                high_boundary = num

        nodes[high_boundary].rule.boundary = level

        self.__set_boundary(nodes,level+1,start,high_boundary+1)
        self.__set_boundary(nodes,level+1,high_boundary,end)

    def get_result(self):
        return self.nodes

    def write_file(self,filename):
        root = et.Element('GPR')
        part = et.SubElement(root,'part')
        part.set('id','P1')
        group = et.SubElement(part,'group')
        self.__write_gpr(group,self.nodes,1)

        document = md.parseString(et.tostring(root,'utf-8'))

        file = open(filename,'w')
        document.writexml(file, encoding='utf-8', newl='\n', indent='', addindent='  ')
        file.close()

    def __write_gpr(self,root,nodes,n):
        divide_flag = 0
        divide_num = 0
        group = et.SubElement(root,'group')
        for i in range(len(nodes)):
            if nodes[i].rule.boundary == n:
                divide_flag = 1
                divide_num = i
            
        if divide_flag == 1:
            self.__write_gpr(group,nodes[:divide_num+1],n+1)
            self.__write_rule(nodes[divide_num],group)
            self.__write_gpr(group,nodes[divide_num:],n+1)
        else:
            for i in range(len(nodes)-1):
                if nodes[i] != None:
                    if nodes[i+1] != None:
                        note = et.SubElement(group,'note')
                        note.set('id',nodes[i].rule.id)
                        if i+1 != len(nodes)-1:
                            self.__write_rule(nodes[i+1],group)

    def __write_rule(self,node,group):
        rules_list = ["2a","2b","3c","3b","3c","3d","4","5","6"]
        for k in rules_list:
            if node.rule.GPR[k] == 1:
                applied = et.SubElement(group,'applied')
                applied.set('rule',k)

    def min_max(self,x):
        ary = np.array(x)
        min = ary.min()
        max = ary.max()
        result = (ary-min)/(max-min)
        return result
        

    
        


    
    
