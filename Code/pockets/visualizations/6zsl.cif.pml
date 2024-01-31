
        from pymol import cmd,stored
        
        set depth_cue, 1
        set fog_start, 0.4
        
        set_color b_col, [36,36,85]
        set_color t_col, [10,10,10]
        set bg_rgb_bottom, b_col
        set bg_rgb_top, t_col      
        set bg_gradient
        
        set  spec_power  =  200
        set  spec_refl   =  0
        
        load "data/6zsl.cif", protein
        create ligands, protein and organic
        select xlig, protein and organic
        delete xlig
        
        hide everything, all
        
        color white, elem c
        color bluewhite, protein
        #show_as cartoon, protein
        show surface, protein
        #set transparency, 0.15
        
        show sticks, ligands
        set stick_color, magenta
        
        
        
        
        # SAS points
 
        load "data/6zsl.cif_points.pdb.gz", points
        hide nonbonded, points
        show nb_spheres, points
        set sphere_scale, 0.2, points
        cmd.spectrum("b", "green_red", selection="points", minimum=0, maximum=0.7)
        
        
        stored.list=[]
        cmd.iterate("(resn STP)","stored.list.append(resi)")    # read info about residues STP
        lastSTP=stored.list[-1] # get the index of the last residue
        hide lines, resn STP
        
        cmd.select("rest", "resn STP and resi 0")
        
        for my_index in range(1,int(lastSTP)+1): cmd.select("pocket"+str(my_index), "resn STP and resi "+str(my_index))
        for my_index in range(1,int(lastSTP)+1): cmd.show("spheres","pocket"+str(my_index))
        for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_scale","0.4","pocket"+str(my_index))
        for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_transparency","0.1","pocket"+str(my_index))
        
        
        
        set_color pcol1 = [0.361,0.576,0.902]
select surf_pocket1, protein and id [7208,6748,6539,6540,6541,6545,6546,6511,6517,6518,6519,7751,6529,6547,6551,6553,6781,6552,6779,6528,7747,6520,7738,7214,7215,7216,8485,7210,7231,7415,7217,7213,6726,6729,7207,7446,7444,7445,8487,8488,8489,8724,7750,6722,6723,6749,8504,8501,8502,8503,8480,8486,6724,6725] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.349,0.702]
select surf_pocket2, protein and id [2334,2354,4060,4067,4069,3022,4066,4304,2116,2122,3327,2331,2151,2123,3326,4083,2385,2386,4081,4082,4062,2796,2797,2798,2799,2789,2790,2791,2150,2133,2144,2145,2991,2141,2134,2137,2792,2795] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.404,0.361,0.902]
select surf_pocket3, protein and id [3034,3041,3044,3047,3036,3037,4043,3030,3114,4223,4247,4192,4220,4221,4222,1325,1339,1340,4042,4044,4199,4200,4202,4195,4201,4203,1324,1327,3920,3897,3908,4027,4028,4029,4198,3910,3663,3050,3054,3110,3083,3111,3109,3118,3115] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.416,0.278,0.702]
select surf_pocket4, protein and id [7468,7471,7472,5796,8620,8622,8615,8619,8086,8087,8082,8324,8326,5767,5785,5798,5773,8461,8462,8463,8464,8667,5797,8313,8623,8447,8448,8336,8329] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.663,0.361,0.902]
select surf_pocket5, protein and id [7529,7534,7535,7538,8643,7527,7465,7474,7479,7507,7539,7542,7461,7590,8656,8658,7557,7558,8640,8641,8642,8657,7460,7468,8620,8612,8614,8615,8667] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.616,0.278,0.702]
select surf_pocket6, protein and id [3042,3046,2693,1356,1357,1342,1353,1336,1081,1046,1047,1048,1043,1044,3066,3067,1056,1082,1083,1084,1682,1348,1350,1027,1025,1026,2832,2833,2860,2861,1024] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.361,0.878]
select surf_pocket7, protein and id [5461,6129,6135,4529,5428,5399,5401,5402,5460,5426,5432,5433,4525,4528,4677,4679,4522,4655,4644,4664,4667,4541,4540,5438] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.278,0.584]
select surf_pocket8, protein and id [7249,7250,5476,7252,5446,7256,7259,5445,5357,5358,5333,5366,7571,7583,7481,7482,7260,5418,5448,5449,5452,7264] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.361,0.620]
select surf_pocket9, protein and id [1705,1720,1721,1722,1723,2850,2852,2882,1727,1738,2902,2881,2689,2687,2690,2910,2911,2909] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.278,0.380]
select surf_pocket10, protein and id [5172,5175,5180,5182,5166,5168,5170,4592,4836,4829,4843,4844,4860,4862,5017,5015,4986,4999,5181] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.361,0.361]
select surf_pocket11, protein and id [5513,5514,5515,5542,5501,5540,7490,7491,5484,5485,7250,7251,7470,5502,5504,5505,5506,7111,7278,7466,7279] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.702,0.380,0.278]
select surf_pocket12, protein and id [321,324,507,490,675,486,506,474,328,316,347,349,663,664] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.620,0.361]
select surf_pocket13, protein and id [2831,2832,2842,2846,2834,2838,3058,3159,3057,3147,2841,991,994,1018,875,900,871,905,899] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.584,0.278]
select surf_pocket14, protein and id [3314,2385,2386,2408,2151,2158,2156,2157,2159,2160,2161,2134,1943] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.902,0.878,0.361]
select surf_pocket15, protein and id [3809,3810,4488,3776,3831,3549,3550,3551,3533] 
set surface_color,  pcol15, surf_pocket15 
   
        
        deselect
        
        orient
        