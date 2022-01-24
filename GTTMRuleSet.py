from abc import ABCMeta, abstractmethod
import numpy as np

class GTTMRuleSet(metaclass=ABCMeta):
    def set_score(self,score):
        self.score = score

    def apply_rules(self):
        nodes = self.get_node()
        self.__apply_rules(self.nodes)

    def __apply_rules(self,nodes):
        if len(nodes) <= 1:
            return 

        rules = self.create_rule(nodes)
        param = self.get_param()

        apply_rule_list = self.rule_set(rules)
        for key,rule_func in apply_rule_list.item():
            for i in range(len(nodes)):
                nodes[i].rule[key] = rule_func(i)
        
        D_sum = np.array([0.0] * len(nodes))

        for i in range(len(nodes)):
            for key in apply_rule_list.keys():
                D_sum[i] += nodes[i].rule[key] * param["S"+key] 
        D_sum[i] += self.apply_recursive_process()

        next_node = self.calc_analysis(D_sum)
        self.__apply_rules(next_node)
                
    
    @abstractmethod
    def get_param(self):
        pass

    @abstractmethod
    def create_rule(self,nodes):
        pass

    @abstractmethod
    def get_nodes(self):
        pass

    @abstractmethod
    def calc_analysis(self):
        pass     

    @abstractmethod
    def rule_set(self):
        pass

    def apply_recursive_process(self):
        return 0

    def get_result(self):
        return self.root

    @abstractmethod
    def write_file(self, filename):
        pass
