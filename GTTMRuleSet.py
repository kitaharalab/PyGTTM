from abc import ABCMeta, abstractmethod
import xml.etree.ElementTree as et
import xml.dom.minidom as md
import numpy as np

class GTTMRuleSet(metaclass=ABCMeta):
    def set_score(self,score):
        self.score = score

    def apply_rules(self):
        nodes = self.get_nodes()
        self.__apply_rules(nodes)

    def __apply_rules(self,nodes):
        if len(nodes) <= 1:
            return 
        
        param = self.get_param()
        rules = self.create_rule(nodes,param)

        apply_rule_list = self.rule_set(rules)
        for key,rule_func in apply_rule_list.items():
            for i in range(len(nodes)):
                nodes[i].rule[key] = rule_func(i)
        
        D_sum = np.array([0.0] * len(nodes))

        for i in range(len(nodes)):
            for key in apply_rule_list.keys():
                D_sum[i] += nodes[i].rule[key] * param["S"+key] 
        D_sum = self.apply_recursive_process(D_sum)

        next_node = self.calc_analysis(nodes,D_sum)
        self.recursion_process(next_node)
    
    def recursion_process(self, nodes):
        self.__apply_rules(nodes)

    @abstractmethod
    def get_param(self):
        pass

    @abstractmethod
    def create_rule(self,nodes,param):
        pass

    def get_nodes(self):
        return self.nodes

    @abstractmethod
    def calc_analysis(self,nodes,D_sum):
        pass     

    @abstractmethod
    def rule_set(self,rules):
        pass

    def apply_recursive_process(self,D_sum):
        return D_sum

    def get_result(self):
        return self.get_nodes()

    def write_file(self, filename):
        root,element = self.set_element()

        self.construct_xml(element)

        document = md.parseString(et.tostring(root,'utf-8'))
        file = open(filename,'w')
        document.writexml(file, encoding='utf-8', newl='\n', indent='', addindent='  ')
        file.close()

    @abstractmethod
    def set_element(self):
        pass

    @abstractmethod
    def construct_xml(self,element):
        pass
