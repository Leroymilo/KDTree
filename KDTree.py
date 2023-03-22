#This script is an implementation of a Kd-tree, used to quickly search the point among a scatter plot closer to a given position .
#Reference : https://youtu.be/BK5x7IUTIyU

#A lot of things would be nicer if I had the dimension each point divides as an attribute of the Node it's in : no need to pass it down the recursive call.
#I might change the code to do that some day but I don't really have time for now, just ideas.

from __future__ import annotations
from typing import Iterable

try :
    import matplotlib.pyplot as plt
    can_plot = True
except ImportError as e :
    print("Warning : matplotlib module not found, draw functions won't work !")
    can_plot = False

class Point :
    def __init__(self, coordinates:tuple[float], dimensions:int=2) :
        #You can add stuff here, this is the minimum needed.
        self.C = coordinates
        self.d = dimensions
    
    def __str__(self) :
        return str(self.C)
    
    def __eq__(self, __o: Point) -> bool:
        return self.C == __o.C
    
    def draw(self, col:str='b', mark:str='s', label:str="") :
        # Requires matplotlib
        if can_plot :
            plt.scatter(self.C[0], self.C[1], c=col, marker=mark, label=label)
    
    def dist(self, __o: Point) :
        # Euclidian distance
        return sum([(self.C[i]-__o.C[i])**2 for i in range(self.d)])**0.5


class Scatter :
    def __init__(self, name:str="Scatter", dimensions:int=2) :
        self.name = name
        self.d = dimensions
        self.points = set()
        self.kdt = KdTree()

    def addPoint(self, point:Point) :
        if point.d == self.d :
            self.points.add(point)
        else :
            print("error : point is in the wrong dimension")
    
    def buildKDT(self) :
        if self.kdt.root is not None :
            self.kdt.clean()
        self.kdt.build(self)

    def draw(self) :
        fig = plt.figure(self.name)
        for point in self.points :
            plt.scatter(point.C[0], point.C[1], c='b', marker='s')
    
    def closest(self, point:tuple[float]) :
        mindist = float('inf')
        point = Point(point, len(point))
        closePoint = None
        
        for p in self.points :
            dist = point.dist(p)
            if dist < mindist :
                closePoint = p
                mindist = dist
        return closePoint
    
    def closest_kdt(self, coords:tuple[float]) :
        #First call of a recursive function
        if self.kdt.root is not None :
            self.kdt.root.clean()
        start = self.kdt.root
        closestNode = self.kdt.root
        point = Point(coords, len(coords))
        closeNode, min_dist = self.search_down(point, start, 0, closestNode, point.dist(start.info))
        return closeNode, min_dist

    def search_down(self, point:Point, deeperNode: Node, dim, closeNode, min_dist) :
        deeperNode.searched = True

        #Going down the tree :
        while deeperNode.hasChildren() :
            prev = deeperNode
            if point.C[dim] <= deeperNode.info.C[dim] and deeperNode.lowChild is not None :
                deeperNode = deeperNode.lowChild
            else :
                deeperNode = deeperNode.highChild

            dist = point.dist(deeperNode.info)
            deeperNode.searched = True
            if dist < min_dist :
                min_dist = dist
                closeNode = deeperNode
            
            dim = (dim+1)%self.d
        
        closeNode, min_dist = self.search_up(point, deeperNode, dim, closeNode, min_dist)
        return closeNode, min_dist
    
    def search_up(self, point:Point, deeperNode: Node, dim, closeNode, min_dist):
        #Going up the tree :
        toCheck: Node = deeperNode.p
        dim = (dim-1)%self.d
        while toCheck is not None :
            if abs(point.C[dim]-toCheck.info.C[dim]) < min_dist :
                #Check on the other side of the line if it's closer than the current minDist
                if toCheck.lowChild is not None and not toCheck.lowChild.searched :
                    closeNode, min_dist = self.search_down(point, toCheck.lowChild, (dim+1)%self.d, closeNode, min_dist)
                elif toCheck.highChild is not None and not toCheck.highChild.searched :
                    closeNode, min_dist = self.search_down(point, toCheck.highChild, (dim+1)%self.d, closeNode, min_dist)

            #Continue going up the tree
            toCheck = toCheck.p
            dim = (dim-1)%self.d
        return closeNode, min_dist


class Node :
    def __init__(self, point:Point, id_:int) :
        self.info = point
        self.id = id_

        self.clear()
    
    def __str__(self) :
        return f'{self.info} [{self.lowChild}, {self.highChild}]'
    
    def __eq__(self, __o: Node) -> bool:
        return self.info == __o.info

    def clear(self) :
        self.p: Node = None
        self.p_side: bool = None

        self.lowChild: Node = None
        self.highChild: Node = None

        self.dim_cut = None

        self.searched = False

        self.balance = 0
        self.height = 1

    def set_child(self, side:bool, child:Node) :
        "side == 0 for lowChild, side == 1 for highChild"
        if side :
            self.highChild = child
        else :
            self.lowChild = child
        
        child.p = self
        child.p_side = side

        self.calc_balance()

    def calc_balance(self) :
        old_h, old_b = self.height, self.balance

        low_h = 0 if self.lowChild is None else self.lowChild.height
        high_h = 0 if self.highChild is None else self.highChild.height

        self.height = max(low_h, high_h) + 1
        self.balance = high_h - low_h

        if self.p is None :
            return None
        if old_b != self.balance or old_h != self.height :
            return self.p.calc_balance()
        
    def build_KDTree(self, lower:list[Node], higher:list[Node], dim:int) :
        self.dim_cut = dim
        nextDim = (dim+1)%self.info.d

        if len(lower) != 0 :
            low_lower, low_med, low_higher = quick_select(lower, (len(lower)-1)//2, nextDim)
            low_med.build_KDTree(low_lower, low_higher, nextDim)

        if len(higher) != 0 :
            high_lower, high_med, high_higher = quick_select(higher, (len(higher)-1)//2, nextDim)
            high_med.build_KDTree(high_lower, high_higher, nextDim)
        
        self.set_child(side=0, child=low_med)
        self.set_child(side=1, child=high_med)
    
    def hasChildren(self) :
        return (self.highChild is not None or self.lowChild is not None)


class KdTree :
    def __init__(self) :
        self.root = None
        self.nodes = []
        self.N = 0

    def __str__(self):
        return str(self.root)

    def build(self, scatter:Iterable[Point]) :
        #Build a balanced tree from an iterable of points, as described in the video linked line 2.
        self.nodes = [Node(point, i) for i, point in enumerate(scatter)]
        self.N = len(self.nodes)
        lower, median, higher = quick_select(self.nodes, (self.N-1)//2, 0)
        self.root = median
        self.root.dim_cut = 0

        self.root.build_KDTree(lower, higher, 1%scatter.d)
    
    def insert(self, point: Point, keep_balance=True) :
        self.nodes.append(Node(point, self.N))
        self.N += 1

        # Inserting a new point in the KD tree as described here : 

#Recursive function used to build the KD tree.
def quick_select(nodes:list[Node], i, dim) -> tuple[list[Node], Node, list[Node]] :
    #Find which point is in the ith position on the dimension dim.
    #Returns the list on lower points, the searched point and the list of higher points.
    pivot: Node = nodes[i]
    lower = []
    higher = []
    for point in nodes[:i] + nodes[i+1:] :
        if point.info.C[dim] <= pivot.info.C[dim] :
            lower.append(point)
        else :
            higher.append(point)
    
    if len(lower) < i :
        reclow, ith, rechigh = quick_select(higher, i-len(lower)-1, dim)
        return lower+[pivot]+reclow, ith, rechigh
    
    elif len(lower) > i :
        reclow, ith, rechigh = quick_select(lower, i, dim)
        return reclow, ith, rechigh+[pivot]+higher
    
    else :
        return lower, pivot, higher