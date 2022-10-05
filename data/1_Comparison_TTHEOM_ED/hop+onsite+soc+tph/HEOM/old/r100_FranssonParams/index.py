import numpy as np
import itertools
import sys

class Hierarchy_index():
    def __init__(self, number_el_state, number_spin, number_lead, number_pole, head_tail_coup, truncation_tier):
        """
        Generating a index list indices[],indices_plus[],indices_minus[] 
        """
        self.number_el_state = number_el_state
        self.number_spin = number_spin
        self.number_lead = number_lead
        self.number_pole = number_pole
        self.head_tail_coup = head_tail_coup
        self.truncation_tier = truncation_tier
        self.modes_imagtime()

    def modes_imagtime(self):
        """
        Generating the index_mode: mode in {[0:number_el_state,0:number_lead,0:number_pole+1,0:sign_pm, sin/cos]}
        i.e., sin/cos=(0,1)[s*cos+(1-s)*sin], number_lead=2(left:0, right:1),    number_el_state=1,  number_pole=1,  sign_pm=2(0:+, 1:-)
        """
        el_states=np.arange(self.number_el_state)
        spins=np.arange(self.number_spin)
        leads=np.arange(self.number_lead)
        poles=np.arange(self.number_pole)
        sign=np.arange(2)
        sin_cos=np.arange(2) #[0,1] 0:sin, 1:cos
    
        modes_imag = []
        if self.head_tail_coup:
            for elem in itertools.product(spins, leads, poles, sign, sin_cos):
                modes_imag.append(list(elem))
        else:
            for elem in itertools.product(el_states, spins, leads, poles, sign, sin_cos):
                modes_imag.append(list(elem))
        self.number_modes_imag = len(modes_imag)
    
        return np.array(modes_imag)

    def indices_real_generation(self):
        self.number_mode_real = self.number_modes_imag // 2
        self.indices_real = []
        self.index_plus_real = []
        tier_range = [0]

        i_tier=0
        self.indices_real.append(np.zeros(self.number_mode_real,dtype=bool))
        tier_range.append(1)

        i_tier = 1
        i_count = 1
        sublist_plus = np.full((self.number_mode_real),-1,dtype=int)
        for i_mode in range(self.number_mode_real):
            index_temp=np.zeros(self.number_mode_real,dtype=bool)
            index_temp[i_mode]=True
            self.indices_real.append(index_temp)
            sublist_plus[i_mode]=i_count
            i_count += 1
        self.index_plus_real.append(sublist_plus)
        tier_range.append(i_count)

        while ( i_tier < self.truncation_tier):
            i_tier += 1
            for i_index in range(tier_range[i_tier-1], tier_range[i_tier]):
                sublist_plus = np.full((self.number_mode_real),-1,dtype=int)
                self.index_plus_real.append(sublist_plus)
            for i_index in range(tier_range[i_tier-1], tier_range[i_tier]):
                index_where=np.where(self.indices_real[i_index])[0]
                for i_mode in range(self.number_mode_real):
                    if i_mode not in index_where:
                        index_temp=self.indices_real[i_index].copy()
                        index_temp[i_mode]=True
                        positions=[i_mode]+list(index_where)
                        positions.sort()
                        hierarchy_plus=0
                        for i_pos in positions:
                            hierarchy_plus=self.index_plus_real[hierarchy_plus][i_pos]
                        if (hierarchy_plus == -1):
                            self.indices_real.append( index_temp )
                            self.index_plus_real[i_index][i_mode]=i_count
                            i_count += 1
                        else:
                            self.index_plus_real[i_index][i_mode]=hierarchy_plus
            tier_range.append(i_count)

        return np.array(self.indices_real)

    def indices_imag_generation(self):
        self.indices_imag = []
        self.index_plus_imag = []
        tier_range = [0]


        i_tier=0
        self.indices_imag.append(np.zeros(self.number_modes_imag,dtype=bool))
        tier_range.append(1)

        i_tier = 1
        i_count = 1
        sublist_plus = np.full((self.number_modes_imag),[-1],dtype=int)
        for i_mode in range(self.number_modes_imag):
            index_temp=np.zeros(self.number_modes_imag,dtype=bool)
            index_temp[i_mode]=True
            self.indices_imag.append(index_temp)
            sublist_plus[i_mode]=i_count
            i_count += 1
        self.index_plus_imag.append(sublist_plus)
        tier_range.append(i_count)

        while ( i_tier < self.truncation_tier):
            i_tier += 1
            for i_index in range(tier_range[i_tier-1], tier_range[i_tier]):
                sublist_plus = np.full((self.number_modes_imag),[-1],dtype=int)
                self.index_plus_imag.append(sublist_plus)
            for i_index in range(tier_range[i_tier-1], tier_range[i_tier]):
                index_where=np.where(self.indices_imag[i_index])[0]
                for i_mode in range(self.number_modes_imag):
                    if i_mode not in index_where:
                        index_temp=self.indices_imag[i_index].copy()
                        index_temp[i_mode]=True
                        positions=[i_mode]+list(index_where)
                        positions.sort()
                        hierarchy_plus=0
                        for i_pos in positions:
                            hierarchy_plus=self.index_plus_imag[hierarchy_plus][i_pos]
                        if (hierarchy_plus == -1):
                            self.indices_imag.append( index_temp )
                            self.index_plus_imag[i_index][i_mode]=i_count
                            i_count += 1
                        else:
                            self.index_plus_imag[i_index][i_mode]=hierarchy_plus
            tier_range.append(i_count)

        return np.array(self.indices_imag)

if __name__=="__main__":
    Number_lead=2
    Number_pole=2 
    Number_spin = 1
    Number_el_state=1

    hierarchy=Hierarchy_index(Number_el_state, Number_spin, Number_lead,Number_pole,2)
    modes = hierarchy.modes_imagtime()
    output_modes=open('modes.txt','w')
    for i in range(len(modes)):
        output_modes.write('index_modes['+str(i)+']\t=\t'+str(modes[i])+'\n')
    output_modes.close()        

    indices=hierarchy.indices_real_generation()
    output_indices=open('indices.txt','w')
    for i in range(len(indices)):
        output_indices.write('indices['+str(i)+']\t= '+str(indices[i])+'\n')
    output_indices.close()

    indices=hierarchy.indices_imag_generation()
    output_indices=open('indices_imag.txt','w')
    for i in range(len(indices)):
        output_indices.write('indices['+str(i)+']\t= '+str(indices[i])+'\n')
    output_indices.close()
