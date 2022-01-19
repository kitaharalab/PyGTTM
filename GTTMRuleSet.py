from abc import ABCMeta, abstractmethod

class GTTMRuleSet(metaclass=ABCMeta):
    def set_score(self,score):
        self.score = score

    @abstractmethod
    def apply_rules(self):
        pass

    @abstractmethod
    def get_result(self):
        pass

    @abstractmethod
    def write_file(self, filename):
        pass
