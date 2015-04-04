n = [];
a = [];
b = [];
c = [];
q_mag = [];
selected_process = [];
q_c_l_x = [];
q_c_l_y = [];
q_c_l_z = [];
q_s_l_x = [];
q_s_l_y = [];
q_s_l_z = [];
q_c_v_x = [];
q_c_v_y = [];
q_c_v_z = [];
q_s_v_x = [];
q_s_v_y = [];
q_s_v_z = [];
theta_q = [];
phi_q = [];
Gamma_0 = [];
self_Scatter_Rate = [];
dGamma_Longitudinal = [];
dGamma_Transverse = [];

selected_process = [];
timeStep = [];

globalTime = [];

valley_i = [];
valley_f = [];

k_beforeField_c_l_x = [];
k_beforeField_c_l_y = [];
k_beforeField_c_l_z = [];
k_beforeField_s_l_x = [];
k_beforeField_s_l_y = [];
k_beforeField_s_l_z = [];
k_beforeField_c_v_x = [];
k_beforeField_c_v_y = [];
k_beforeField_c_v_z = [];
k_beforeField_s_v_x = [];
k_beforeField_s_v_y = [];
k_beforeField_s_v_z = [];

grp_v_beforeField_c_l_x = [];
grp_v_beforeField_c_l_y = [];
grp_v_beforeField_c_l_z = [];
grp_v_beforeField_s_l_x = [];
grp_v_beforeField_s_l_y = [];
grp_v_beforeField_s_l_z = [];
grp_v_beforeField_c_v_x = [];
grp_v_beforeField_c_v_y = [];
grp_v_beforeField_c_v_z = [];
grp_v_beforeField_s_v_x = [];
grp_v_beforeField_s_v_y = [];
grp_v_beforeField_s_v_z = [];

pos_beforeField_c_l_x = [];
pos_beforeField_c_l_y = [];
pos_beforeField_c_l_z = [];
pos_beforeField_s_l_x = [];
pos_beforeField_s_l_y = [];
pos_beforeField_s_l_z = [];
pos_beforeField_c_v_x = [];
pos_beforeField_c_v_y = [];
pos_beforeField_c_v_z = [];
pos_beforeField_s_v_x = [];
pos_beforeField_s_v_y = [];
pos_beforeField_s_v_z = [];

kineticEnergy_beforeField = [];

k_afterField_c_l_x = [];
k_afterField_c_l_y = [];
k_afterField_c_l_z = [];
k_afterField_s_l_x = [];
k_afterField_s_l_y = [];
k_afterField_s_l_z = [];
k_afterField_c_v_x = [];
k_afterField_c_v_y = [];
k_afterField_c_v_z = [];
k_afterField_s_v_x = [];
k_afterField_s_v_y = [];
k_afterField_s_v_z = [];

grp_v_afterField_c_l_x = [];
grp_v_afterField_c_l_y = [];
grp_v_afterField_c_l_z = [];
grp_v_afterField_s_l_x = [];
grp_v_afterField_s_l_y = [];
grp_v_afterField_s_l_z = [];
grp_v_afterField_c_v_x = [];
grp_v_afterField_c_v_y = [];
grp_v_afterField_c_v_z = [];
grp_v_afterField_s_v_x = [];
grp_v_afterField_s_v_y = [];
grp_v_afterField_s_v_z = [];

pos_afterField_c_l_x = [];
pos_afterField_c_l_y = [];
pos_afterField_c_l_z = [];
pos_afterField_s_l_x = [];
pos_afterField_s_l_y = [];
pos_afterField_s_l_z = [];
pos_afterField_c_v_x = [];
pos_afterField_c_v_y = [];
pos_afterField_c_v_z = [];
pos_afterField_s_v_x = [];
pos_afterField_s_v_y = [];
pos_afterField_s_v_z = [];

kineticEnergy_afterField = [];

k_afterScatter_c_l_x = [];
k_afterScatter_c_l_y = [];
k_afterScatter_c_l_z = [];
k_afterScatter_s_l_x = [];
k_afterScatter_s_l_y = [];
k_afterScatter_s_l_z = [];
k_afterScatter_c_v_x = [];
k_afterScatter_c_v_y = [];
k_afterScatter_c_v_z = [];
k_afterScatter_s_v_x = [];
k_afterScatter_s_v_y = [];
k_afterScatter_s_v_z = [];

grp_v_afterScatter_c_l_x = [];
grp_v_afterScatter_c_l_y = [];
grp_v_afterScatter_c_l_z = [];
grp_v_afterScatter_s_l_x = [];
grp_v_afterScatter_s_l_y = [];
grp_v_afterScatter_s_l_z = [];
grp_v_afterScatter_c_v_x = [];
grp_v_afterScatter_c_v_y = [];
grp_v_afterScatter_c_v_z = [];
grp_v_afterScatter_s_v_x = [];
grp_v_afterScatter_s_v_y = [];
grp_v_afterScatter_s_v_z = [];

pos_afterScatter_c_l_x = [];
pos_afterScatter_c_l_y = [];
pos_afterScatter_c_l_z = [];
pos_afterScatter_s_l_x = [];
pos_afterScatter_s_l_y = [];
pos_afterScatter_s_l_z = [];
pos_afterScatter_c_v_x = [];
pos_afterScatter_c_v_y = [];
pos_afterScatter_c_v_z = [];
pos_afterScatter_s_v_x = [];
pos_afterScatter_s_v_y = [];
pos_afterScatter_s_v_z = [];

kineticEnergy_afterScatter = [];



fid = fopen('Data.txt');
tline = fgetl(fid);


for n = 1:11
    tline = fgetl(fid);
end

tline = fgetl(fid);
disp(tline)

while ~strcmp(tline, 'END DATA')
    %fgetl(fid);
    fgetl(fid);
    fgetl(fid);

    n = [n str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    a = [a str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    b = [b str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    c = [c str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);
    
    q_mag = [q_mag str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    selected_process = [selected_process str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    q_c_l_x = [q_c_l_x str2double(fgetl(fid))];
    q_c_l_y = [q_c_l_y str2double(fgetl(fid))];
    q_c_l_z = [q_c_l_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    q_s_l_x = [q_s_l_x str2double(fgetl(fid))];
    q_s_l_y = [q_s_l_y str2double(fgetl(fid))];
    q_s_l_z = [q_s_l_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    q_c_v_x = [q_c_v_x str2double(fgetl(fid))];
    q_c_v_y = [q_c_v_y str2double(fgetl(fid))];
    q_c_v_z = [q_c_v_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    q_s_v_x = [q_s_v_x str2double(fgetl(fid))];
    q_s_v_y = [q_s_v_y str2double(fgetl(fid))];
    q_s_v_z = [q_s_v_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    theta_q = [theta_q str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);
    
    phi_q = [phi_q str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    Gamma_0 = [Gamma_0 str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    self_Scatter_Rate = [self_Scatter_Rate str2double(fgetl(fid))];
    fgetl(fid);    
    fgetl(fid);

    dGamma_Longitudinal = [dGamma_Longitudinal str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);
    
    dGamma_Transverse = [dGamma_Transverse str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);

    selected_process = [selected_process str2double(fgetl(fid))];
    fgetl(fid);

    timeStep = [timeStep str2double(fgetl(fid))];    
    fgetl(fid);
    fgetl(fid);

    globalTime = [globalTime str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    valley_i = [valley_i str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    valley_f = [valley_f str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    k_beforeField_c_l_x = [k_beforeField_c_l_x str2double(fgetl(fid))];
    k_beforeField_c_l_y = [k_beforeField_c_l_y str2double(fgetl(fid))];
    k_beforeField_c_l_z = [k_beforeField_c_l_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    k_beforeField_s_l_x = [k_beforeField_s_l_x str2double(fgetl(fid))];
    k_beforeField_s_l_y = [k_beforeField_s_l_y str2double(fgetl(fid))];
    k_beforeField_s_l_z = [k_beforeField_s_l_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    k_beforeField_c_v_x = [k_beforeField_c_v_x str2double(fgetl(fid))];
    k_beforeField_c_v_y = [k_beforeField_c_v_y str2double(fgetl(fid))];
    k_beforeField_c_v_z = [k_beforeField_c_v_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    k_beforeField_s_v_x = [k_beforeField_s_v_x str2double(fgetl(fid))];
    k_beforeField_s_v_y = [k_beforeField_s_v_y str2double(fgetl(fid))];
    k_beforeField_s_v_z = [k_beforeField_s_v_z str2double(fgetl(fid))];
    
    fgetl(fid);
    fgetl(fid);
    
    grp_v_beforeField_c_l_x = [grp_v_beforeField_c_l_x str2double(fgetl(fid))];
    grp_v_beforeField_c_l_y = [grp_v_beforeField_c_l_y str2double(fgetl(fid))];
    grp_v_beforeField_c_l_z = [grp_v_beforeField_c_l_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    grp_v_beforeField_s_l_x = [grp_v_beforeField_s_l_x str2double(fgetl(fid))];
    grp_v_beforeField_s_l_y = [grp_v_beforeField_s_l_y str2double(fgetl(fid))];
    grp_v_beforeField_s_l_z = [grp_v_beforeField_s_l_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    grp_v_beforeField_c_v_x = [grp_v_beforeField_c_v_x str2double(fgetl(fid))];
    grp_v_beforeField_c_v_y = [grp_v_beforeField_c_v_y str2double(fgetl(fid))];
    grp_v_beforeField_c_v_z = [grp_v_beforeField_c_v_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    grp_v_beforeField_s_v_x = [grp_v_beforeField_s_v_x str2double(fgetl(fid))];
    grp_v_beforeField_s_v_y = [grp_v_beforeField_s_v_y str2double(fgetl(fid))];
    grp_v_beforeField_s_v_z = [grp_v_beforeField_s_v_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);
    
    pos_beforeField_c_l_x = [pos_beforeField_c_l_x str2double(fgetl(fid))];
    pos_beforeField_c_l_y = [pos_beforeField_c_l_y str2double(fgetl(fid))];
    pos_beforeField_c_l_z = [pos_beforeField_c_l_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    pos_beforeField_s_l_x = [pos_beforeField_s_l_x str2double(fgetl(fid))];
    pos_beforeField_s_l_y = [pos_beforeField_s_l_y str2double(fgetl(fid))];
    pos_beforeField_s_l_z = [pos_beforeField_s_l_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    pos_beforeField_c_v_x = [pos_beforeField_c_v_x str2double(fgetl(fid))];
    pos_beforeField_c_v_y = [pos_beforeField_c_v_y str2double(fgetl(fid))];
    pos_beforeField_c_v_z = [pos_beforeField_c_v_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    pos_beforeField_s_v_x = [pos_beforeField_s_v_x str2double(fgetl(fid))];
    pos_beforeField_s_v_y = [pos_beforeField_s_v_y str2double(fgetl(fid))];
    pos_beforeField_s_v_z = [pos_beforeField_s_v_z str2double(fgetl(fid))];
    fgetl(fid);
        
    kineticEnergy_beforeField = [kineticEnergy_beforeField str2double(fgetl(fid))];

    fgetl(fid);
    fgetl(fid);
    
    k_afterField_c_l_x = [k_afterField_c_l_x str2double(fgetl(fid))];
    k_afterField_c_l_y = [k_afterField_c_l_y str2double(fgetl(fid))];
    k_afterField_c_l_z = [k_afterField_c_l_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    k_afterField_s_l_x = [k_afterField_s_l_x str2double(fgetl(fid))];
    k_afterField_s_l_y = [k_afterField_s_l_y str2double(fgetl(fid))];
    k_afterField_s_l_z = [k_afterField_s_l_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    k_afterField_c_v_x = [k_afterField_c_v_x str2double(fgetl(fid))];
    k_afterField_c_v_y = [k_afterField_c_v_y str2double(fgetl(fid))];
    k_afterField_c_v_z = [k_afterField_c_v_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    k_afterField_s_v_x = [k_afterField_s_v_x str2double(fgetl(fid))];
    k_afterField_s_v_y = [k_afterField_s_v_y str2double(fgetl(fid))];
    k_afterField_s_v_z = [k_afterField_s_v_z str2double(fgetl(fid))];

    fgetl(fid);
    fgetl(fid);
    
    grp_v_afterField_c_l_x = [grp_v_afterField_c_l_x str2double(fgetl(fid))];
    grp_v_afterField_c_l_y = [grp_v_afterField_c_l_y str2double(fgetl(fid))];
    grp_v_afterField_c_l_z = [grp_v_afterField_c_l_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    grp_v_afterField_s_l_x = [grp_v_afterField_s_l_x str2double(fgetl(fid))];
    grp_v_afterField_s_l_y = [grp_v_afterField_s_l_y str2double(fgetl(fid))];
    grp_v_afterField_s_l_z = [grp_v_afterField_s_l_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    grp_v_afterField_c_v_x = [grp_v_afterField_c_v_x str2double(fgetl(fid))];
    grp_v_afterField_c_v_y = [grp_v_afterField_c_v_y str2double(fgetl(fid))];
    grp_v_afterField_c_v_z = [grp_v_afterField_c_v_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    grp_v_afterField_s_v_x = [grp_v_afterField_s_v_x str2double(fgetl(fid))];
    grp_v_afterField_s_v_y = [grp_v_afterField_s_v_y str2double(fgetl(fid))];
    grp_v_afterField_s_v_z = [grp_v_afterField_s_v_z str2double(fgetl(fid))];

    fgetl(fid);
    fgetl(fid);
    
    pos_afterField_c_l_x = [pos_afterField_c_l_x str2double(fgetl(fid))];
    pos_afterField_c_l_y = [pos_afterField_c_l_y str2double(fgetl(fid))];
    pos_afterField_c_l_z = [pos_afterField_c_l_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);


    pos_afterField_s_l_x = [pos_afterField_s_l_x str2double(fgetl(fid))];
    pos_afterField_s_l_y = [pos_afterField_s_l_y str2double(fgetl(fid))];
    pos_afterField_s_l_z = [pos_afterField_s_l_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    pos_afterField_c_v_x = [pos_afterField_c_v_x str2double(fgetl(fid))];
    pos_afterField_c_v_y = [pos_afterField_c_v_y str2double(fgetl(fid))];
    pos_afterField_c_v_z = [pos_afterField_c_v_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    pos_afterField_s_v_x = [pos_afterField_s_v_x str2double(fgetl(fid))];
    pos_afterField_s_v_y = [pos_afterField_s_v_y str2double(fgetl(fid))];
    pos_afterField_s_v_z = [pos_afterField_s_v_z str2double(fgetl(fid))];

    fgetl(fid);
    
    kineticEnergy_afterField = [kineticEnergy_afterField str2double(fgetl(fid))];
    
    fgetl(fid);
    fgetl(fid);
    
    k_afterScatter_c_l_x = [k_afterScatter_c_l_x str2double(fgetl(fid))];
    k_afterScatter_c_l_y = [k_afterScatter_c_l_y str2double(fgetl(fid))];
    k_afterScatter_c_l_z = [k_afterScatter_c_l_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    k_afterScatter_s_l_x = [k_afterScatter_s_l_x str2double(fgetl(fid))];
    k_afterScatter_s_l_y = [k_afterScatter_s_l_y str2double(fgetl(fid))];
    k_afterScatter_s_l_z = [k_afterScatter_s_l_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    k_afterScatter_c_v_x = [k_afterScatter_c_v_x str2double(fgetl(fid))];
    k_afterScatter_c_v_y = [k_afterScatter_c_v_y str2double(fgetl(fid))];
    k_afterScatter_c_v_z = [k_afterScatter_c_v_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    k_afterScatter_s_v_x = [k_afterScatter_s_v_x str2double(fgetl(fid))];
    k_afterScatter_s_v_y = [k_afterScatter_s_v_y str2double(fgetl(fid))];
    k_afterScatter_s_v_z = [k_afterScatter_s_v_z str2double(fgetl(fid))];

    fgetl(fid);
    fgetl(fid);

    grp_v_afterScatter_c_l_x = [grp_v_afterScatter_c_l_x str2double(fgetl(fid))];
    grp_v_afterScatter_c_l_y = [grp_v_afterScatter_c_l_y str2double(fgetl(fid))];
    grp_v_afterScatter_c_l_z = [grp_v_afterScatter_c_l_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    grp_v_afterScatter_s_l_x = [grp_v_afterScatter_s_l_x str2double(fgetl(fid))];
    grp_v_afterScatter_s_l_y = [grp_v_afterScatter_s_l_y str2double(fgetl(fid))];
    grp_v_afterScatter_s_l_z = [grp_v_afterScatter_s_l_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    grp_v_afterScatter_c_v_x = [grp_v_afterScatter_c_v_x str2double(fgetl(fid))];
    grp_v_afterScatter_c_v_y = [grp_v_afterScatter_c_v_y str2double(fgetl(fid))];
    grp_v_afterScatter_c_v_z = [grp_v_afterScatter_c_v_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    grp_v_afterScatter_s_v_x = [grp_v_afterScatter_s_v_x str2double(fgetl(fid))];
    grp_v_afterScatter_s_v_y = [grp_v_afterScatter_s_v_y str2double(fgetl(fid))];
    grp_v_afterScatter_s_v_z = [grp_v_afterScatter_s_v_z str2double(fgetl(fid))];

    fgetl(fid);
    fgetl(fid);
    
    pos_afterScatter_c_l_x = [pos_afterScatter_c_l_x str2double(fgetl(fid))];
    pos_afterScatter_c_l_y = [pos_afterScatter_c_l_y str2double(fgetl(fid))];
    pos_afterScatter_c_l_z = [pos_afterScatter_c_l_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);


    pos_afterScatter_s_l_x = [pos_afterScatter_s_l_x str2double(fgetl(fid))];
    pos_afterScatter_s_l_y = [pos_afterScatter_s_l_y str2double(fgetl(fid))];
    pos_afterScatter_s_l_z = [pos_afterScatter_s_l_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    pos_afterScatter_c_v_x = [pos_afterScatter_c_v_x str2double(fgetl(fid))];
    pos_afterScatter_c_v_y = [pos_afterScatter_c_v_y str2double(fgetl(fid))];
    pos_afterScatter_c_v_z = [pos_afterScatter_c_v_z str2double(fgetl(fid))];
    fgetl(fid);
    fgetl(fid);

    pos_afterScatter_s_v_x = [pos_afterScatter_s_v_x str2double(fgetl(fid))];
    pos_afterScatter_s_v_y = [pos_afterScatter_s_v_y str2double(fgetl(fid))];
    pos_afterScatter_s_v_z = [pos_afterScatter_s_v_z str2double(fgetl(fid))];

    fgetl(fid);
    
    kineticEnergy_afterScatter = [kineticEnergy_afterScatter str2double(fgetl(fid))];
    fgetl(fid);
    
    tline = fgetl(fid);
    disp(tline)
end

fclose(fid);
