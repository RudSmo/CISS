module GlobalConstants

# PHYSICAL CONSTANTS [kg m s A Cd]

global const hbar() = 1.0545718*10^(-34)  
global const elma() = 9.10938356*10^(-31) 
global const gyro() = 176085964400.0
global const timc() = 2.41888*10^(-17)


# MISC OBJECTS

global ex = [1,0,0]
global ey = [0,1,0]
global ez = [0,0,1]
global s_x = [0. 1; 1 0]
global s_y = [0. -1im; 1im 0]
global s_z = [1. 0; 0 -1]

global AppendData = 100
global SaveData = Int(AppendData * 10)

end
