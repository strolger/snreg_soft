define pts
expand 0.25
points
expand 1
end

define doit
data xform.out
xcolumn 3
ycolumn 4
ucolumn 5
vcolumn 6
limits
box
pts
set u u - x
set v v - y
vfield 1 100
end
