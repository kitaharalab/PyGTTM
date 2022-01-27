import numpy as np
import xml.etree.ElementTree as et
import xml.dom.minidom as md
from GTTMRuleSet import GTTMRuleSet

class MPR_node:

    def __init__(self, L_end, velo, valu, vol, slur, num, id):
        self.L_end = L_end #拍点の位置
        self.velo = velo #拍点から始まる音のベロシティ
        self.valu = valu #音価
        self.vol = vol #連続する音量の長さ
        self.slur = slur #スラーの長さ
        self.num = num #音高
        self.boundary = 0 #start:1
        self.id = id
        self.dot = 1
        self.rule = {"1":0, "2":0, "3":0, "4":0, "5a":0, "5b":0, "5c":0, "5d":0, "5e":0}


class Rules:
    def __init__(self,nodes,velo,valu,vol,slur,num):
        self.nodes = nodes
        self.velo = velo
        self.valu = valu
        self.vol = vol
        self.slur = slur
        self.num = num

    def MPR1(self, i, k):
        i_start,i_end = self.group_set(i)
        x = self.MPR1_x(i,k,i_start,i_end)
        y = self.MPR1_y(i,k,i_start,i_end)
        z = self.MPR1_z(i,k,i_start,i_end)

        if x == 0 or y == 0:
            return 0
        if (y/x * self.Wr) + (z/y * (1-self.Wr)) > self.T1:
            return 1
        else :
            return 0

    def MPR1_x(self, i, k, i_start, i_end):
        i_sum = 0
        for j in range(i_start, i_end):
            i_sum += 1 if self.velo[j] > 0 else 0

        k_start,k_end = self.group_set(k)
        k_sum = 0
        for l in range(k_start, k_end):
            k_sum += 1 if self.velo[l] > 0 else 0

        return i_sum + k_sum

    def MPR1_y(self, i, k, i_start, i_end):
        sum = 0
        for j in range(i_start, i_end):
            sum += 1 if (j-i+k) < len(self.nodes) and self.velo[i]>0 and self.velo[j-i+k]>0 else 0

        return sum

    def MPR1_z(self, i, k, i_start, i_end):
        sum = 0
        for j in range(i_start, i_end):
            sum += 1 if (j-i+k) < len(self.nodes) and self.velo[i] > 0 and self.num[j-1] == self.num[j] and self.num[k+j-i-1] == self.num[k+j-i] else 0

        return sum

    def group_set (self, i):
        k = i
        while self.nodes[k].boundary == 0:
            k -= 1
        
        if i < len(self.nodes)-1:
            j = i+1
        else:
            j = i
        while self.nodes[j].boundary == 0 and j != len(self.nodes)-1:
            j += 1
        
        return k,j

       #最も強い拍がグループの中で比較的早く出る拍節構造を優先する.
    def MPR2(self, i):
        i_start,i_end = self.group_set(i)
        #print(i,i_start,i_end)
        if i_end-1 == 0 or i_end-i_start == 0:
            return 0
        return (i_end - i)/(i_end - i_start)

    #拍点に音符がある拍節構造を優先する.
    def MPR3(self, i):
        if self.velo[i] > 0:
            return 1
        else:
            return 0

    #強く弾いた拍が強拍となる拍節構造を優先する
    def MPR4(self, i):
        if self.velo[i] > 2 * np.average(self.velo) * self.T4:
            return 1
        else :
            return 0

    #5a,5b,5cでは相対的に長い音，長い音量，長いスラーが強拍となる拍節構造を優先する.
    def MPR5a(self, i):
        if self.valu[i] > 2 * np.average(self.valu) * self.T5a:
            return 1
        else :
            return 0

    def MPR5b(self, i):
        if self.vol[i] > 2 * np.average(self.vol) * self.T5b:
            return 1
        else :
            return 0

    def MPR5c(self, i):
        if self.slur[i] > 2 * np.average(self.slur) * self.T5c:
            return 1
        else :
            return 0

    #相対的に長い アーティキュレーションパターンの繰り返し が強拍となる拍節構造を優先する.
    def MPR5d(self, i):
        if i < len(self.nodes)-1 and self.nodes[i].rule["5a"] == 1 and self.nodes[i+1].rule["5a"] == 1:
            return 1
        else :
            return 0

    #同一音高が連続している場合に強拍となる拍節構造を優先する．
    def MPR5e(self, i):
        if i < len(self.nodes)-1 and self.num[i] == self.num[i+1]:
            return 1
        else :
            return 0

    def set_param(self,param):
        self.Wr = param["Wr"]
        self.T1 = param["T1"]
        self.T4 = param["T4"]
        self.T5a = param["T5a"]
        self.T5b = param["T5b"]
        self.T5c = param["T5c"]
        

class MPR(GTTMRuleSet):
    param = {"Wr":0.5,"S1":0.5,"S2":0.5,"S3":0.5,"S4":0.5,"S5a":0.5,"S5b":0.5,"S5c":0.5,"S5d":0.5,
            "S5e":0.5,"S10":0.5,"T1":0.5,"T4":0.5,"T5a":0.5,"T5b":0.5,"T5c":0.5}

    def __init__(self,score,GPR,**kwargs):
        self.set_score(score)
        self.nodes = self.__set_MPR(GPR)
        for key,value in kwargs:
            self.param[key] = value

    def __set_MPR(self, GPR):
        k = 0
        M_nodes = []
        velo = 0
        valu = 0
        vol = 0
        slur = 0
        num = 0
        div = self.score.div
        node_len = self.score.node_len
        for i in range(node_len):
            L_end = i / div
            if k >= self.score.get_score_length() or self.score.get_L_end(k) != L_end :
                M_nodes.append(MPR_node(L_end, 0, 0, 0, 0, 0, ""))
            else:
                velo = self.score.get_velo(k)
                valu = self.score.get_valu(k)
                vol = self.score.get_vol(k,L_end)
                slur = self.score.get_slur(k)
                num = self.score.get_num(k)
                id = self.score.get_id(k)
                node = MPR_node(L_end, velo, valu, vol, slur, num, id)
                
                if GPR[k].boundary != 0:
                    node.boundary = 1
                M_nodes.append(node)
                k += 1
        
        for n in M_nodes:
            if n.velo == 1:
                vol = n.vol
                num = n.num
            else:
                vol -= 1/self.div
                if vol > 0.0:
                    n.vol = vol
                    n.num = num
        
        M_nodes[0].boundary = 1
        return M_nodes


    def get_param(self):
        return self.param

    def create_rule(self,nodes,param):
        velo = []
        valu = []
        vol = []
        slur = []
        num = []
        for i in range(len(nodes)):
            velo.append(nodes[i].velo)
            valu.append(nodes[i].valu)
            vol.append(nodes[i].vol)
            slur.append(nodes[i].slur)
            num.append(nodes[i].num)

        self.rules = Rules(nodes,velo,valu,vol,slur,num)
        self.rules.set_param(self.param)
        return self.rules

    def rule_set(self,rules):
        apply_rule = {"2":rules.MPR2,"3":rules.MPR3,"4":rules.MPR4,"5a":rules.MPR5a,"5b":rules.MPR5b,"5c":rules.MPR5c,"5e":rules.MPR5e,"5d":rules.MPR5d}
        return apply_rule

    def apply_recursive_process(self,D_sum):
        D_low = np.array([0.0] * len(D_sum))
        for i in range(len(D_sum)):
            sum = 0
            for k in range(len(D_sum)):
                if self.rules.MPR1(i,k) == 1:
                    sum += D_sum[k] * self.param["S1"]
            D_low[i] += D_sum[i] + sum
        return D_low

    def calc_analysis(self,nodes,D_sum):
        next_m = [0,0,0,0,0]
        for i in range(len(nodes)):
            for m in range(5):
                if (m==0 or m==1) and np.mod(i-m,2) == 0:
                    next_m[m] += D_sum[i]
                if (m==2 or m==3 or m==4) and np.mod(i-m,3) == 2:
                    next_m[m] += D_sum[i] * self.param["S10"]

        m = np.argmax(next_m)
        next_nodes = []
        for i in range(len(nodes)):
            if (m == 0 or m == 1) and np.mod(i-m,2) == 0:
                nodes[i].dot += 1
                if i >= 1 and nodes[i-1].boundary == 1:
                    nodes[i].boundary = 1
                if len(next_nodes) > 0 and next_nodes[-1].boundary == 1:
                    nodes[i].boundary = 0
                next_nodes.append(nodes[i])

            if (m == 2 or m == 3 or m == 4) and np.mod(i-m,3) == 0:
                nodes[i].dot += 1
                if (i >= 1 and nodes[i-1].boundary == 1) or (i >= 2 and nodes[i-2].boundary == 1):
                    nodes[i].boundary = 1
                if len(next_nodes) > 0 and next_nodes[-1].boundary == 1:
                    nodes[i].boundary = 0
                next_nodes.append(nodes[i])

        return next_nodes


    def get_result(self):
        return self.nodes

    def set_element(self):
        root = et.Element('MPR')
        part = et.SubElement(root,'part')
        part.set('id','P1')
        return root, part

    def construct_xml(self, element):
        for n in self.nodes:
            metric = et.SubElement(element,'metric')
            metric.set('dot',str(n.rule.dot))
            metric.set('at',str(n.L_end))
            if n.rule.id != "":
                note = et.SubElement(metric,'note')
                note.set('id', n.rule.id)

