function new_world  = addObject( world, object, insert_pos)
%
%   function addObject.m
%
%   Martin Krucinski
%   Adds point cloud objec to point cloud world
%   requires both point clouds to be at the same pitch_x and pitch_y
%   but not necessarily same size
%
%   2019-03-13  Inital version
%               object_insert_pos is an x and y index (2D vector)
%               (not a position [m] yet)

size_world  = size(world.Z);
new_world   = world;

size_object = size(object.Z);
NWX         = size_world(2);
NWY         = size_world(1);

NOX         = size_object(2);
NOY         = size_object(1);

for x=1:NOX,
    for y=1:NOY,
        insert_x   = x + insert_pos(1);
        insert_y   = y + insert_pos(2);
        
        if (insert_x<=NWX) && (insert_y<=NWY),
            new_world.Z(insert_y, insert_x)   = max( new_world.Z(insert_y, insert_x) , object.Z(y,x) );
        end
    end
end

        
        
