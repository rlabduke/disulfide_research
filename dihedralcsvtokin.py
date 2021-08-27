#1B80,B:CYS::   3: ,B:CYS::  15: , 301.91, 298.84, 279.06, 300.04, 295.76
#1B80,B:CYS::  14: ,B:CYS:: 285: , 300.09, 295.02, 282.83, 305.24, 292.64
#1B80,B:CYS::  34: ,B:CYS:: 120: , 179.51, 104.85, 85.46, 282.83, 306.60
"""
@kinemage
@group{group name}
@dotlist{list name}
"""
csvfile = open("updateddihedrals.csv", "r")
print("@kinemage")
print("@dimensions {chi1} {chi2} {chi3} {chi2prime} {chi1prime}")
print("@dimminmax 0 360 0 360 0 360 0 360 0 360")
print("@group{group name} dimension=5 select")
print("@dotlist{list name} dimension=5")
for line in csvfile: 
    x = line.split(",")
    pointid = "{" + x[0] + x[1] + x[2] + "}"
    coordinates = " ".join(x[3:])
    
    print(pointid + coordinates.rstrip())