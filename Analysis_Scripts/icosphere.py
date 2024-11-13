import numpy as np
from scipy.spatial.transform import Rotation as R

def vertex(x, y, z): 
    """ Return vertex coordinates fixed to the unit sphere """ 
    length = np.sqrt(x**2 + y**2 + z**2) 
    return [i / length for i in (x,y,z)] 

def middle_point(verts,middle_point_cache,point_1, point_2): 
    """ Find a middle point and project to the unit sphere """ 
    # We check if we have already cut this edge first 
    # to avoid duplicated verts 
    smaller_index = min(point_1, point_2) 
    greater_index = max(point_1, point_2) 
    key = '{0}-{1}'.format(smaller_index, greater_index) 
    if key in middle_point_cache: return middle_point_cache[key] 
    # If it's not in cache, then we can cut it 
    vert_1 = verts[point_1] 
    vert_2 = verts[point_2] 
    middle = [sum(i)/2 for i in zip(vert_1, vert_2)] 
    verts.append(vertex(*middle)) 
    index = len(verts) - 1 
    middle_point_cache[key] = index 
    return index

def icosphere(subdiv,radius):
    # verts for icosahedron
    r = (1.0 + np.sqrt(5.0)) / 2.0;
    verts = np.array([[-1.0, r, 0.0],[ 1.0, r, 0.0],[-1.0, -r, 0.0],
                      [1.0, -r, 0.0],[0.0, -1.0, r],[0.0, 1.0, r],
                      [0.0, -1.0, -r],[0.0, 1.0, -r],[r, 0.0, -1.0],
                      [r, 0.0, 1.0],[ -r, 0.0, -1.0],[-r, 0.0, 1.0]]);
    # rescale the size to radius of 0.5
    verts /= np.linalg.norm(verts[0])
    # adjust the orientation
    r = R.from_quat([[0.19322862,-0.68019314,-0.19322862,0.68019314]])
    verts = r.apply(verts)
    verts = list(verts)
    # indicates vertices that form each of the triangles
    faces = [[0, 11, 5],[0, 5, 1],[0, 1, 7],[0, 7, 10],
             [0, 10, 11],[1, 5, 9],[5, 11, 4],[11, 10, 2],
             [10, 7, 6],[7, 1, 8],[3, 9, 4],[3, 4, 2],
             [3, 2, 6],[3, 6, 8],[3, 8, 9],[5, 4, 9],
             [2, 4, 11],[6, 2, 10],[8, 6, 7],[9, 8, 1],];
    
    for i in range(subdiv):
        middle_point_cache = {}
        faces_subdiv = []
        for tri in faces: 
            v1  = middle_point(verts,middle_point_cache,tri[0], tri[1])
            v2  = middle_point(verts,middle_point_cache,tri[1], tri[2])
            v3  = middle_point(verts,middle_point_cache,tri[2], tri[0])
            faces_subdiv.append([tri[0], v1, v3]) 
            faces_subdiv.append([tri[1], v2, v1]) 
            faces_subdiv.append([tri[2], v3, v2]) 
            faces_subdiv.append([v1, v2, v3]) 
        faces = faces_subdiv
                
    return [np.array(verts)*radius/1.0, np.array(faces)]

def generateunitnormals(meshvar,radius):
    vertices=meshvar[0];
    triangles=meshvar[1];
    tsize=np.shape(triangles)
    unitnormal=[]
    frac=radius*0.9
    for i in range(0,tsize[0]):
        unitnormtriangle=[];
        face = triangles[i,:]
        #get the vertex coordinates
        v1 = vertices[face[0],:]
        v2 = vertices[face[1],:]
        v3 = vertices[face[2],:]
        # generate vectors
        vec1 = v2 - v1
        vec2 = v3 - v2
        vec3 = v1 - v3
        vec4 = v1 - frac*v1
        vec5 = v2 - frac*v2
        vec6 = v3 - frac*v3
        # +ve distance is outside
        #cross product
        cpvec12=np.cross(vec1,vec2)
        normval=np.linalg.norm(cpvec12)
        cpvec12=cpvec12/normval
        d = cpvec12[0]*v1[0] + cpvec12[1]*v1[1] + cpvec12[2]*v1[2]
        cpvec12 = np.append(cpvec12,d)

        cpvec14=np.cross(vec1,vec4)
        normval=np.linalg.norm(cpvec14)
        cpvec14=cpvec14/normval
        d = cpvec14[0]*v1[0] + cpvec14[1]*v1[1] + cpvec14[2]*v1[2]
        cpvec14 = np.append(cpvec14,d)

        cpvec25=np.cross(vec2,vec5)
        normval=np.linalg.norm(cpvec25)
        cpvec25=cpvec25/normval
        d = cpvec25[0]*v2[0] + cpvec25[1]*v2[1] + cpvec25[2]*v2[2]
        cpvec25 = np.append(cpvec25,d)

        cpvec36=np.cross(vec3,vec6)
        normval=np.linalg.norm(cpvec36)
        cpvec36=cpvec36/normval
        d = cpvec36[0]*v1[0] + cpvec36[1]*v1[1] + cpvec36[2]*v1[2]
        cpvec36 = np.append(cpvec36,d)

        unitnormtriangle.append(cpvec12)
        unitnormtriangle.append(cpvec14)
        unitnormtriangle.append(cpvec25)
        unitnormtriangle.append(cpvec36)
        unitnormal.append(np.array(unitnormtriangle))
    return unitnormal
# https://computergraphics.stackexchange.com/questions/10831/how-can-i-get-a-signed-distance-sdf-from-a-mesh

# Create bead triangle neighbor list with a simple filter
# Key bead , value - triangle
# Use normal to find distance
# If a triangle is within range, check the neighboring triangles.
# Assign the bead to just one neighbor
# Keep a counter in each triangle
def generateNeighborList(meshvar):
    # Simpler implementation
    #Step 1.
    # Find neighbors of each triangle into a list of lists
    # Go through the vertices
    # For each vertex, find triangles that contain the vertex.
    vertices = meshvar[0]
    triangles = meshvar[1]
    Nvert = (np.shape(vertices))[0]
    Ntriangles = (np.shape(triangles))[0]
    VTneighborlist = [] #Vertex triangle neighborlist
    TTneighborlist= [[]] * Ntriangles
    for i in range(0,Nvert):
        vtosearch  = i
        neighbors = []
        for j in range(0,Ntriangles):
            tvec = triangles[j][:]
            if(np.in1d(vtosearch,tvec)):
                neighbors.append(j)
        neighbors = np.sort(neighbors)
        VTneighborlist.append(neighbors)
        # Go through the neighbors. Add them to each other's list
        for t1iter, t1 in enumerate(neighbors):   
            temp = TTneighborlist[t1]     
            temp2 = np.setdiff1d(neighbors,t1).tolist()
            if(len(temp)):            
                temp2 = temp+temp2
                temp2.sort()
            TTneighborlist[t1] = temp2
    return [TTneighborlist,VTneighborlist]

def interpolateallmonomers (fc):
    interpcoord = []
    Nbeads = (np.shape(fc))[0]
    interpcoords = []
    for i in range(0,Nbeads-1):
        cyl1 = fc[i,:]
        cyl2 = fc[i+1,:]
        cyllength = np.linalg.norm(cyl2-cyl1)
        nmonomers = int(np.ceil(cyllength/2.7e-3))
        #Add minus end
        interpcoords.append(cyl1)
        # Add rest of the monomers
        for mid in range(1,nmonomers-1):
            alpha = mid/nmonomers
            tmp = cyl1*(1-alpha) + alpha*cyl2
            interpcoords.append(tmp)
        #Add plus end
        interpcoords.append(cyl2)
    return np.array(interpcoords)
        
#Go through triangle
#Step 2.
# Assign a triangle for each bead
# For first bead in the filament, go through all triangles
# For other beads, start with triangle that prev. bead is assigned to, then its neighbors before going on to other triangles.
from enum import Enum
class SearchAlgoType(Enum):
    ORIGINAL = 1
    LISTEDSEARCH = 2
    SPHERICALSEARCH = 3

# Check 2A - distance from top triangle plane
def checkdist_bead_tangle(bc,t,donestatus, unitnormal,actincounter, printStatus=False):
    if(donestatus[t]):
        return False
    else:
        donestatus[t] = True
    foundstatus = False
    dthreshold = 0.2
    normmat = unitnormal[t] 
    unvec = normmat[0,:] #Has 4 normals

    distval = abs(unvec[0]*bc[0] + unvec[1]*bc[1] + unvec[2]*bc[2] - unvec[3])
    if(distval<dthreshold):
        # Check 2B - If the distance constraint is satisfied,
        #check if point is within the cake slice

        slice1norm = normmat[1,:]
        slice2norm = normmat[2,:]
        slice3norm = normmat[3,:]
        slice1dist = 0
        slice2dist = 0
        slice3dist = 0
        for d in range(0,3):
            slice1dist = slice1dist + slice1norm[d]*bc[d]
            slice2dist = slice2dist + slice2norm[d]*bc[d]
            slice3dist = slice3dist + slice3norm[d]*bc[d]
        if(slice1dist <=0 and slice2dist <=0 and slice3dist <=0):
            if(printStatus):
                print('Distance from center='+str(np.linalg.norm(bc)))
                print('Distance from triangle plane='+str(distval))
            foundstatus = True
            actincounter[t] = actincounter[t] + 1;
    return foundstatus

def generatedensityfield(meshvar, unitnormal, filcoord, searchType = SearchAlgoType.ORIGINAL, Nlists = [], sphericalc=[]):
    printStatus = False
    vertices = meshvar[0]
    triangles = meshvar[1]
    Nvert = (np.shape(vertices))[0]
    Ntriangles = (np.shape(triangles))[0]
    dthreshold = 0.2
    dthreshold_4m_center =  0.9*0.9
    actincounter = list(np.zeros((Ntriangles,)))
    total_actin  = 0
    for f, fc in enumerate(filcoord):
        # Interpolate each monomer in the filament
        interpcoord = interpolateallmonomers(fc)
        #interpcoord = fc
        Nbeads = (np.shape(interpcoord))[0]
        counter = 0
        init_tangle_list = []

        if searchType==SearchAlgoType.SPHERICALSEARCH:
            tanthetavec = sphericalc[0]
            cosphivec = sphericalc[1]

        for b in range(0,Nbeads):
            foundstatus = False
            total_actin = total_actin + 1
            bc = interpcoord[b,:]
            # Check 1 - distance from center
            dfromcenter = bc[0]*bc[0]+bc[1]*bc[1]+bc[2]*bc[2]
            if(dfromcenter>=dthreshold_4m_center):
                #tracks if this triangle has been considered for this particular bead
                donestatus = [False] * Ntriangles
                
                t_search_list = list(np.arange(0,Ntriangles))

                if searchType==SearchAlgoType.LISTEDSEARCH:
                    TTNlist = Nlists[0];
                    #go through initial list
                    if init_tangle_list:
                        for t in init_tangle_list:
                            # Check 2A - distance from top triangle plane
                            foundstatus=checkdist_bead_tangle(bc,t,donestatus, unitnormal,actincounter);
                            if(foundstatus):
                                init_tangle_list = [t]+TTNlist[t]
                                break;
                    if not(foundstatus):
                        for tier, t in enumerate(t_search_list):
                            foundstatus=checkdist_bead_tangle(bc,t,donestatus, unitnormal,actincounter);
                            if(foundstatus):
                                init_tangle_list = [t]+TTNlist[t]
                                break;
                                
                elif searchType==SearchAlgoType.SPHERICALSEARCH:
                    costhetavec = sphericalc[0]
                    cosphivec = sphericalc[1]
                    Rgtangle = 0.027
                    zthreshold = 1/Rgtangle

                else:
                    for t in range(0,Ntriangles):
                        if(donestatus[t]):
                            continue;
                        else:
                            donestatus[t] = 1

                        # Check 2A - distance from top triangle plane
                        normmat = unitnormal[t] 
                        unvec = normmat[0,:] #Has 4 normals

                        distval = abs(unvec[0]*bc[0] + unvec[1]*bc[1] + unvec[2]*bc[2] - unvec[3])

                        if(distval<dthreshold):
                            # Check 2B - If the distance constraint is satisfied,
                            #check if point is within the cake slice

                            slice1norm = normmat[1,:]
                            slice2norm = normmat[2,:]
                            slice3norm = normmat[3,:]
                            slice1dist = 0
                            slice2dist = 0
                            slice3dist = 0
                            for d in range(0,3):
                                slice1dist = slice1dist + slice1norm[d]*bc[d]
                                slice2dist = slice2dist + slice2norm[d]*bc[d]
                                slice3dist = slice3dist + slice3norm[d]*bc[d]
                            if(slice1dist <=0 and slice2dist <=0 and slice3dist <=0):
                                if(printStatus):
                                    foundtriangle=t
                                    print('Found '+str([f,b]))
                                    print('Distance from center='+str(np.linalg.norm(bc)))
                                    print('Distance from triangle plane='+str(distval))
                                foundstatus = True
                                actincounter[t] = actincounter[t] + 1;
                                break;
                    
            if(not(foundstatus) and printStatus):
                print('Could not find '+str([f,b]))
                print('Distance from center='+str(np.linalg.norm(bc)))
    return np.array(actincounter)