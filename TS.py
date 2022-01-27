from os import register_at_fork
import numpy as np
import xml.etree.ElementTree as et
import xml.dom.minidom as md
from GTTMRuleSet import GTTMRuleSet

class TS_node:
    def __init__(self, L_end, R_end, rest, ioi, pitch, dot, boundary, id):
        self.head_L_end = L_end #ヘッドの位置
        self.head_R_end = R_end
        self.ts_len = R_end-L_end
        self.L_end = L_end
        self.R_end = R_end
        self.rest = rest #次のヘッドまでの消音時間から発音時間の差
        self.ioi = ioi #次のヘッドまでの発音時間から発音時間までの差
        self.pitch = pitch #次のヘッドとの音高差
        self.dot = dot #現在のヘッドの拍点
        self.boundary = boundary #境界のレベル
        self.id = id
        self.primary = None
        self.secondary = None
        self.left = None
        self.right = None
        self.prebrunch = None
        self.d = 0
        self.rule = {"1":0, "3":0, "4":0, "8":0, "9":0 }

class TS_rule:

    def __init__(self):
        self.TS = {"1":0, "3":0, "4":0, "8":0, "9":0 }
        self.id = ""

class Rules:
    def __init__(self,nodes,rest,ioi,pitch,dot):
        self.nodes = nodes
        self.rest = rest
        self.ioi = ioi
        self.pitch = pitch 
        self.dot = dot
    #拍点の数が多いものを選好する.
    def TS1(self,i):
        return self.dot[i] / np.max(self.dot)

    #音高が高い音，バスの低い音を優先する
    def TS3(self,i):
        if np.min(self.pitch) == np.max(self.pitch):
            return 0
        else:
            high_pitch = (self.pitch[i] - np.min(self.pitch))/(np.max(self.pitch) - np.min(self.pitch))
            low_pitch = 1- high_pitch

            if high_pitch >= 0.8:
                return high_pitch
            elif low_pitch >= 0.8:
                return low_pitch
            else:
                return 0

    #平行構造を尊重する
    def TS4(self,i,k):
        if i == 0 or i == len(self.nodes)-1 or k == 0 or k == len(self.nodes)-1:
            return 0
        if self.ioi[i-1] == self.ioi[k-1] and self.ioi[i] == self.ioi[k] and self.ioi[i+1] == self.ioi[k+1]:
            return 1
        else :
            return 0
    
    #構造の開始点を選好する
    def TS8(self,i):
        if self.nodes[i].boundary == 1:
            return 1
        else:
            return 0

    #構造の終点を選好する
    def TS9(self,i):
        if i == len(self.nodes)-1 or self.nodes[i+1].boundary == 1:
            return 1
        else:
            return 0


class TS(GTTMRuleSet):
    param = {"Wm":0.5,"Wl":0.5,"Ws":0.5,"S1":0.5,"S3":0.5,"S4":0.5,"S8":0.5,"S9":0.5}
    root = None

    def __init__(self,score,GPR,MPR,**kwargs):
        self.set_score(score)
        self.nodes = self.__set_TS(GPR,MPR)
        for key,value in kwargs:
            self.param[key] = value

    def __set_TS(self,GPR, MPR):
        k = 0
        T_nodes = []
        L_end = 0
        ioi = 0
        rest = 0
        dot = 0
        pitch = 0
        boundary = 0
        id = ""
        boundary_level = 0
        for node in GPR:
            if boundary_level < node.boundary:
                boundary_level = node.boundary
        boundary_level += 1
        for i in range(self.score.get_score_length()):
            L_end = self.score.get_L_end(i)
            R_end = self.score.get_R_end(i)
            if i != self.score.get_score_length()-1:
                rest = self.score.get_rest(i)
                ioi = self.score.get_ioi(i)
            else:
                rest = 0
                ioi = 0
            while MPR[k].L_end < L_end:
                k += 1
            dot = MPR[k].dot
            pitch = self.score.get_num(i)
            
            if GPR[i].boundary != 0:
                if i != 0:
                    boundary = boundary_level - GPR[i].boundary
            else :
                boundary = 0
            id = self.score.get_id(i)
            T_nodes.append(TS_node(L_end,R_end,rest,ioi,pitch,dot,boundary,id))

        T_nodes[0].boundary = boundary_level
        return T_nodes

    def get_param(self):
        return self.param

    def create_rule(self,nodes,param):
        rest = []
        ioi = []
        pitch = []
        dot = []

        for n in nodes:
            rest.append(n.rest)
            ioi.append(n.ioi)
            pitch.append(n.pitch)
            dot.append(n.dot)

        self.rules = Rules(nodes,rest,ioi,pitch,dot)
        return self.rules
        
    def rule_set(self,rules):
        apply_rule = {"1":rules.TS1,"3":rules.TS3,"8":rules.TS8,"9":rules.TS9}
        return apply_rule

    def apply_recursive_process(self, D_sum):
        D_i = np.array([0.0] * len(D_sum))
        for i in range(len(D_sum)):
            sum = 0
            for k in range(len(D_sum)):
                if self.rules.TS4(i,k) == 1:
                    sum += D_sum[i]*self.param["S4"]
            D_i[i] = D_sum[i] + sum
        return D_i

    def calc_analysis(self, nodes, D_sum):
        
        next_brunch = []
        i_end = 0
        i = 0
        tmp_node = None
        D_sum = self.__min_max(D_sum)
        if self.__only_group(nodes) == 1:
            for n in nodes:
                n.boundary -= 1
        #print()
        #for n in nodes:
        #    print(n.head_L_end,n.id,n.boundary)
        while i < len(nodes):
            if i_end == i:
                if tmp_node != None:
                    next_brunch.append(tmp_node)
                    tmp_node = None
                i_start,i_end = self.__group_set(i,nodes)

            #print(i,i_end)
            """
            if D_sum[i] > 0.7:
                if tmp_node != None:
                    next_brunch.append(tmp_node)
                tmp_node = nodes[i]
                i += 1
                print(-9)
            #"""
            
            if i+1 >= i_end :
                if tmp_node != None:
                    if D_sum[i-1] >= D_sum[i]:
                        brunch = TS_node(nodes[i-1].head_L_end,nodes[i-1].head_R_end,nodes[i-1].rest,
                                nodes[i-1].ioi,nodes[i-1].pitch,nodes[i-1].dot,nodes[i-1].boundary,nodes[i-1].id)
                        brunch.primary = nodes[i-1]
                        brunch.secondary = nodes[i]
                        
                    else:
                        brunch = TS_node(nodes[i].head_L_end,nodes[i].head_R_end,nodes[i].rest,
                                nodes[i].ioi,nodes[i].pitch,nodes[i].dot,nodes[i].boundary,nodes[i].id)
                        brunch.primary = nodes[i]
                        brunch.secondary = nodes[i-1]
                        if nodes[i-1].boundary != 0:
                            brunch.boundary = nodes[i-1].boundary
                    if brunch.L_end > brunch.secondary.L_end:
                        brunch.L_end = brunch.secondary.L_end                 
                    if brunch.R_end < brunch.secondary.R_end:
                        brunch.R_end = brunch.secondary.R_end
                    if brunch.L_end > brunch.primary.L_end:
                        brunch.L_end = brunch.primary.L_end
                    if brunch.R_end < brunch.primary.R_end:
                        brunch.R_end = brunch.primary.R_end
                    brunch.ts_len = brunch.R_end-brunch.L_end
                    
                    #print(D_sum[i], D_sum[i+1])
                    #print(brunch.primary.head_L_end,brunch.secondary.head_L_end)
                    next_brunch.append(brunch)
                    tmp_node = None
                    i += 1
                else:
                    next_brunch.append(nodes[i])
                    i += 1
            else:
                if tmp_node != None and (i+1 >= i_end or D_sum[i-1] > D_sum[i+1]):
                    if D_sum[i-1] >= D_sum[i]:
                        brunch = TS_node(nodes[i-1].head_L_end,nodes[i-1].head_R_end,nodes[i-1].rest,
                                nodes[i-1].ioi,nodes[i-1].pitch,nodes[i-1].dot,nodes[i-1].boundary,nodes[i-1].id)
                        brunch.primary = nodes[i-1]
                        brunch.secondary = nodes[i]
                        
                    else:
                        brunch = TS_node(nodes[i].head_L_end,nodes[i].head_R_end,nodes[i].rest,
                                nodes[i].ioi,nodes[i].pitch,nodes[i].dot,nodes[i].boundary,nodes[i].id)
                        brunch.primary = nodes[i]
                        brunch.secondary = nodes[i-1]
                        if nodes[i-1].boundary != 0:
                            brunch.boundary = nodes[i-1].boundary
                    i += 1
                
                else:
                    if D_sum[i] >= D_sum[i+1]:
                        brunch = TS_node(nodes[i].head_L_end,nodes[i].head_R_end,nodes[i].rest,
                                nodes[i].ioi,nodes[i].pitch,nodes[i].dot,nodes[i].boundary,nodes[i].id)
                        brunch.primary = nodes[i]
                        brunch.secondary = nodes[i+1]
                        
                    else:
                        brunch = TS_node(nodes[i+1].head_L_end,nodes[i+1].head_R_end,nodes[i+1].rest,
                                nodes[i+1].ioi,nodes[i+1].pitch,nodes[i+1].dot,nodes[i+1].boundary,nodes[i+1].id)
                        brunch.primary = nodes[i+1]
                        brunch.secondary = nodes[i]
                        if nodes[i].boundary != 0:
                            brunch.boundary = nodes[i].boundary
                    if tmp_node != None:
                        next_brunch.append(tmp_node)
                    i += 2
                
                
                if brunch.L_end > brunch.secondary.L_end:
                    brunch.L_end = brunch.secondary.L_end                 
                if brunch.R_end < brunch.secondary.R_end:
                    brunch.R_end = brunch.secondary.R_end
                if brunch.L_end > brunch.primary.L_end:
                    brunch.L_end = brunch.primary.L_end
                if brunch.R_end < brunch.primary.R_end:
                    brunch.R_end = brunch.primary.R_end
                brunch.ts_len = brunch.R_end-brunch.L_end
                
                next_brunch.append(brunch)
                tmp_node = None

        if tmp_node != None:
            next_brunch.append(tmp_node)
            tmp_node = None

        return next_brunch

    def recursion_process(self, nodes):
        if len(nodes) > 1:
            super().recursion_process(nodes)
        else :
            self.root = nodes[0]


    def __only_group(self, nodes):
        if nodes == None:
            return False
        else:
            is_only = 1
            for n in nodes:
                if n.boundary == 0:
                    is_only = 0
                    break
        return is_only

    def __group_set (self, i,nodes):
        k = i
        while k > 0  and nodes[k].boundary == 0:
            k -= 1
        
        j = i+1
        while j < len(nodes) and nodes[j].boundary == 0 :
            j += 1
        
        return k,j

    def __min_max(self,x):
        ary = np.array(x)
        min = ary.min()
        max = ary.max()
        if min == max:
            return np.zeros(len(x))
        result = (ary-min)/(max-min)
        return result

    def get_result(self):
        return self.root

    def write_file(self,filename):
        root = et.Element('tstree')

        self.__write_ts(root,self.root)

        document = md.parseString(et.tostring(root,'utf-8'))
        file = open(filename,'w')
        document.writexml(file, encoding='utf-8', newl='\n', indent='', addindent='  ')
        file.close()

    def __write_ts(self,root,node):
        if node == None:
            return
        ts = et.SubElement(root,'ts')
        ts.set('timespan',str(node.ts_len))
        ts.set('leftend', str(node.L_end))
        ts.set('rightend',str(node.R_end))
        head = et.SubElement(ts, 'head')
        chord = et.SubElement(head, 'chord')
        chord.set('duration', str(node.head_R_end - node.head_L_end))
        chord.set('velocity', str(90))
        note = et.SubElement(chord, 'note')
        note.set('id',node.id)

        if node.primary != None:
            primary = et.SubElement(ts,'primary')
            self.__write_ts(primary,node.primary)
        if node.secondary != None:
            secondary = et.SubElement(ts, 'secondary')
            self.__write_ts(secondary,node.secondary)
