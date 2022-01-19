import xml.etree.ElementTree as ET
import pretty_midi

class Node:
    def __init__(self, L_end, R_end, a_R_end, Pitch, id):
        self.L_end = L_end
        self.R_end = R_end
        self.act_R_end = a_R_end
        self.pitch = Pitch
        self.velo = 1
        self.id = id

class Score:
    def __init__(self, filename):
        self.filename = filename
        self.nodes = self.__makeScore(filename)
        self.__read_division(filename)

    def get_L_end(self,i):
        return self.nodes[i].L_end
    
    def get_R_end(self,i):
        return self.nodes[i].R_end

    def get_rest(self,i):
        return abs(self.nodes[i+1].L_end - self.nodes[i].R_end)

    def get_ioi(self,i):
        return abs(self.nodes[i+1].L_end - self.nodes[i].L_end)

    def get_regi(self,i):
        return abs(self.nodes[i+1].pitch - self.nodes[i].pitch)

    def get_dyn(self,i):
        #Musicxmlでは読み取れない
        return 0.1

    def get_arti(self,i):
        #return abs((nodes[i+1].R_end-nodes[i+1].L_end)/(nodes[i+1].act_R_end-nodes[i+1].L_end) - (nodes[i].R_end-nodes[i].L_end)/(nodes[i].act_R_end-nodes[i].L_end))
        return 0.1

    def get_leng(self,i):
        return abs((self.nodes[i+1].R_end - self.nodes[i+1].L_end) - (self.nodes[i].R_end - self.nodes[i].L_end))

    def get_id(self,i):
        return self.nodes[i].id

    def get_valu(self,i):
        return self.nodes[i].R_end - self.nodes[i].L_end

    def get_velo(self,i):
        return self.nodes[i].velo

    def get_vol(self,i,L_end):
        tmp = i+1
        if tmp < len(self.nodes):
            while tmp < len(self.nodes) and self.nodes[tmp].L_end == self.nodes[tmp-1].R_end :
                tmp += 1
        return self.nodes[tmp-1].R_end - L_end 

    def get_slur(self,i):
        #スラーの始まりのとき1
        return 0

    def get_num(self,i):
        return self.nodes[i].pitch

    def get_score_length(self):
        return len(self.nodes)

    def __read_division(self, filename):
        score = ET.parse(filename)
        root = score.getroot()
        part = root.find("part")
        self.measure_count = 0
        self.node_len = 0
        for measure in part.findall("measure"):
            self.measure_count += 1
            for attributes in measure.findall("attributes"):
                for divisions in attributes.findall("divisions"):
                    self.div = int(divisions.text)

        for measure in part.findall("measure"):
            for attributes in measure.findall("attributes"):
                for time in attributes.findall("time"):
                    for beats in time.findall("beats"):
                        beat = int(beats.text)
                    for type in time.findall("beat-type"):
                        beat_type = int(type.text)
            self.node_len += 1.0 * beat * (4.0/beat_type) * self.div
        self.node_len = int(self.node_len)

    def __makeScore(self,filename):
        score = ET.parse(filename)
        root = score.getroot()
        part = root.find("part")
        keys = {'C':0,'D':2,'E':4,'F':5,'G':7,'A':9,'B':11}
        nodes = []
        measure_len = 0
        tie_flag = 0
        across_measure = 0
        tmp = 0
        measure_tmp = 0
        note_num = 1
        L_end = 0
        for measure in part.findall("measure"):
            measure_number = int(measure.get("number"))
            if L_end != 0:
                measure_len += L_end
            if measure_number == 0:
                measure_number = 1
            #print("measure :" + str(measure_number))
            for attributes in measure.findall("attributes"):
                for divisions in attributes.findall("divisions"):
                    div = int(divisions.text)
                    #print("div :" + str(div))
                for time in attributes.findall("time"):
                    for beats in time.findall("beats"):
                        beat = int(beats.text)

            if tie_flag != 1:
                L_end = 0.0
            else:
                across_measure = 1
                measure_number -= 1
            key = 0
            for note in measure.findall("note"):
                alt = 0
                for tie in note.findall("tie"):
                    if tie.get("type") == "start":
                        tie_flag = 1
                    elif tie.get("type") == "stop":
                        tie_flag = 2


                for duration in note.findall("duration"):
                    dur = 1.0 / div * int(duration.text)
                    #print(L_end)

                for rest in note.findall("rest"):
                    note_num += 1

                for pitch in note.findall("pitch"):

                    for step in pitch.findall("step"):
                        key = keys[step.text]

                    for alter in pitch.findall("alter"):
                        alt = int(alter.text)

                    for octave in pitch.findall("octave"):
                        key += alt + 12 * int(octave.text)

                    if measure_tmp != measure_number:
                        note_num = 1
                        measure_tmp = measure_number

                    if tie_flag == 0:
                        nodes.append(Node(L_end + measure_len,L_end + measure_len + dur,
                        L_end + measure_len + dur*0.8,key, "P1-"+str(measure_number)+"-"+str(note_num)))
                        note_num += 1
                    elif tie_flag == 2:
                        nodes.append(Node(L_end + measure_len,L_end + measure_len + dur + tmp,
                        L_end + measure_len + dur*0.8 + tmp,key, "P1-"+str(measure_number)+"-"+str(note_num)))
                        note_num += 1
                    

                
                if tie_flag == 0:
                    L_end += dur
                elif tie_flag == 1:
                    tmp += dur
                elif tie_flag == 2:
                    
                    if across_measure == 1:
                        measure_number += 1
                        across_measure = 0
                        L_end = dur
                    else:
                        L_end += dur + tmp 
                    tie_flag = 0
                    tmp = 0

        return nodes