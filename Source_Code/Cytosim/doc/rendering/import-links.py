# This script will import Cytosim's connections into Blender
# It must be run from within Blender, and will read "link####.txt" files
# to generate a collection of cylinders
#
# The input files should be generated by Cytosim's tool "cymart"
# and copied to a local blender directory called 'cymart'
#
# Author: Francois Nedelec, 05.02.2018 -- 06.02.2018

import bpy
from mathutils import Vector
import os.path, math

# Variables
file_path = "cymart/link%04d.txt"
root_name = "link"

scale = 100
deltaZ = 0.030
origin = Vector((0,0,0))
radius = 0.00275 * 0.5 * scale
arp23 = 1

materials = (None, "BioShader-LightRed", "BioShader-DarkRed", "BioShader-Orange", "BioShader-Orange")

# Shortcuts
bld = bpy.data
scene = bpy.context.scene

def remove_objects(name):
    for obj in scene.objects:
        if obj.name.startswith(name):
            bld.objects.remove(obj)


def hide_until(obj, atr, f):
    setattr(obj, atr, True)
    obj.keyframe_insert(atr, frame=f-1)
    setattr(obj, atr, False)
    obj.keyframe_insert(atr, frame=f)


def set_location(obj, pos, f):
    obj.location = pos
    obj.keyframe_insert("location", frame=f)


def set_rotation(obj, rot, f):
    obj.rotation_euler = rot
    obj.keyframe_insert("rotation_euler", frame=f)


def import_frame(frame_id, frame_data):
    """ Import objects from one file"""
    obj_id = 0
    # Go over all the line:
    for line in frame_data:

        data = line.split(" ")
        obj_id = int(data[1])
        obj_type = int(data[2])
        
        pos1 = Vector((float(data[3]), float(data[4]), float(data[5])+deltaZ)) * scale
        pos2 = Vector((float(data[6]), float(data[7]), float(data[8])+deltaZ)) * scale
        
        position = ( pos1 + pos2 ) * 0.5
        length = ( pos1 - pos2 ).length
        diff = ( pos1 - pos2 ) / length       
        ang1 = math.atan2(math.sqrt(diff.x*diff.x+diff.y*diff.y), diff.z);
        ang2 = math.atan2(diff.y, diff.x);
 
        obj = None
        obj_name = root_name + "%04d" % obj_id
        # Find object or create it
        if obj_name in scene.objects:
            # Just get a handle to the existing object
            obj = scene.objects[obj_name]
        else:
            # place cylinder
            mesh = bpy.ops.mesh.primitive_cylinder_add(vertices=16, radius=radius, depth=length)
            obj = bpy.context.object
            obj.name = obj_name
            # place the end for verification
            if 0:
                mesh = bpy.ops.mesh.primitive_cube_add(radius=radius, location=pos1)
                bpy.context.object.name = obj_name+'a'
                mesh = bpy.ops.mesh.primitive_cube_add(radius=radius, location=pos2)
                bpy.context.object.name = obj_name+'b'
            if obj_type == arp23:
                mesh = bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=5, size=2*radius)
                obj2 = bpy.context.object
                obj2.location = Vector((0,0,-length/2.0))
                obj2.name = obj_name+'a'
                obj2.parent = obj
                hide_until(obj2, "hide", frame_id)
                hide_until(obj2, "hide_render", frame_id)
            #scene.objects.link(obj)
            # Hide the object on the previous frames
            hide_until(obj, "hide", frame_id)
            hide_until(obj, "hide_render", frame_id)
        
        # Adjust material:
        if materials[obj_type] in bld.materials:
            obj.active_material = None
            obj.material_slots[obj.active_material_index].link = "OBJECT"
            obj.active_material = bld.materials[materials[obj_type]]

        # set current location
        set_location(obj, position, frame_id)
        set_rotation(obj, (0, ang1, ang2), frame_id)


def import_file(f):
    """ read file number 'f'"""
    file_name = bpy.path.abspath('//') + file_path % f
    # Check if the frame file exists
    if os.path.isfile(file_name):
        #print("reading ", frame_file)
        scene.frame_end = f+1
        # Read the frame file
        with open(file_name, 'r') as file:
            data = file.readlines()
        # Clean line breaks
        data = [line.strip() for line in data]
        import_frame(f+1, data)
    else:
        print("file " + frame_file + " not found")


remove_objects(root_name);

if 0:
    import_file(100)
else:
    # Specify the frames to import here [FIRST, LAST+1]:
    for f in range(0, 101):
        import_file(f)
