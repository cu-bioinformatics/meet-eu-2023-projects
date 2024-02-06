from numpy import prod
import sys
path=sys.argv[1]
size_x=float(sys.argv[2])
size_y=float(sys.argv[3])
size_z=float(sys.argv[4])
search_box=[size_x,size_y,size_z]

method="coordinates"
output="eval"

scale=float(sys.argv[5])

data=open(path, 'r').readlines()
coord = [[float(line[31:38]), float(line[39:46]), float(line[47:54])] for line in data if line.split()[0] in ('ATOM', 'HETATM')]
xcoor, ycoor, zcoor = zip(*coord)
X,Y,Z = [list(xcoor), list(ycoor), list(zcoor)]
ranges=[[min(X), max(X)],[min(Y), max(Y)],[min(Z), max(Z)]]

box = [i*scale for i in [ranges[0][1]-ranges[0][0], ranges[1][1]-ranges[1][0], ranges[2][1]-ranges[2][0]]]

if sum([c==0 for c in box])>1:
    evaluation=False
else:
    evaluation=True
if method in ["volume","vol"]:
    search_box_volume=prod(search_box)
    box_volume=prod(box)
    if box_volume > volume:
        evaluation=False
elif method in ["coord","coordinates"]:
    search_box_sorted=sorted(search_box)
    box_sorted=sorted(box)
    for box_coord, search_box_coord in zip(box_sorted, search_box_sorted):
        if box_coord <= search_box_coord:
            continue
        else:
            evaluation=False
            break

if output in ["evaluation", "eval", "test"]:
    print(evaluation)
elif output in ["size", "molsize"]:
    print(box)