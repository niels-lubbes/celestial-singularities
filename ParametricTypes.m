##############################################################################
#                                                                            #
# Maple 13 computations for the following article by Niels Lubbes            #
#                                                                            #
# Self-intersections of surfaces that contain two circles through each point #
#                                                                            #
#                                                                            #
# To run either copy past the code into a Maple session or execute the       #
# the following Linux command:                                               #
#                                                                            #
# $ maple < ParametricTypes.m                                                #
#                                                                            #
##############################################################################

#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version, see http://www.gnu.org/licenses/
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#

#
# Short description (see articles for details):
#
# A parametric type corresponds to a birational map
# from P1xP1 to a surface X that lies in S3,
# where P1xP1 denotes the product of the projective line with itself and
# S3 denotes the projectivized unit 3-sphere.
# After taking an affine chart, we obtain a map from complex plane C2 to X.
# We compute the components of the preimage in C2 of the singular components in X
# such that a general complex point of this component in X
# has at least two preimages in the parameter domain C2.
#

restart;

kernelopts(printbytes=false);

# global variables
global s,t,u,v,y;

# predefined parametric types
D3a:=[1,-1,0,0,1,0,3/2,0]:
D3b:=[1/10,1/2,0,0,1/2,0,1/2,0]:
D3c:=[3/10,1/2,0,0,1/2,0,1/2,0]:
D3d:=[1,-4,0,0,1,0,0,-1]:
D3e:=[1/5,1/2,0,0,1/2,0,1/2,0]:
D4a:=[1,1,0,0,1,0,1,0]:
D4b:=[1/2,0,1,0,1/2,1,0,0]:
D4c:=[1/2,1/2,0,0,1/2,0,1/2,0]:
D4d:=[1,-3/2,0,0,1,0,-3/2,0]:
D5a:=[1,0,0,0,1,3/2,0,0]:
D5b:=[1,0,0,0,1,2,0,0]:
D5c:=[1,0,0,0,1,5/2,0,0]:

hp := proc(a::list,b::list)
    description "Hamiltonian product of points in the Moebius quadric S3.":
    return [
        a[1]*b[1],
        a[2]*b[2]-a[3]*b[3]-a[4]*b[4]-a[5]*b[5],
        a[3]*b[2]+a[2]*b[3]-a[5]*b[4]+a[4]*b[5],
        a[4]*b[2]+a[5]*b[3]+a[2]*b[4]-a[3]*b[5],
        a[5]*b[2]-a[4]*b[3]+a[3]*b[4]+a[2]*b[5]
    ]
end proc:

isp := proc(p::list)
    description "Inverse stereographic projection P3 --> S3.":
    return [
        p[1]^2+p[2]^2+p[3]^2+p[4]^2,
        2*p[1]*p[2],
        2*p[1]*p[3],
        2*p[1]*p[4],
        -p[1]^2+p[2]^2+p[3]^2+p[4]^2
    ]
end proc:

getPolyList := proc(paramType::list)
    description "Construct ideal for given parametric type paramType.":
    #
    # Example
    #
    # > D4a:=[1,1,0,0,1,0,1,0]:
    # > lprint(getPolyList(D4a));
    #
    local coss,sins,cost,sint,cosu,sinu,cosv,sinv,
          t0,t1,t2,t3,u0,u1,u2,u3,
          mapst,mapuv,
          p1,p2,p3,p4,p5,q1,q2,q3,q4,q5;

    # substitute circle parametrization for cos and sin
    #
    coss, sins := op([(2*s)/(s^2+1),(s^2-1)/(s^2+1)]):
    cost, sint := op(subs(s=t,[coss,sins])):
    cosu, sinu := op(subs(s=u,[coss,sins])):
    cosv, sinv := op(subs(s=v,[coss,sins])):

    # recover map into projective space from parametric type in the vars s,t and u,v
    #
    t0,t1,t2,t3,u0,u1,u2,u3 := op(paramType):
    mapst:=hp(isp([1,t0*coss+t1,t0*sins+t2,t3]),isp([1,u0*cost+u1,u0*sint+u2,u3])):
    mapuv:=hp(isp([1,t0*cosu+t1,t0*sinu+t2,t3]),isp([1,u0*cosv+u1,u0*sinv+u2,u3])):
    mapst:=evala(Simplify(mapst*~((1+s^2)*(1+t^2)))):
    mapuv:=evala(Simplify(mapuv*~((1+u^2)*(1+v^2)))):
    p1,p2,p3,p4,p5:=op(mapst):
    q1,q2,q3,q4,q5:=op(mapuv):

    # Provide ideal that encodes the following set
    #   { (s,t,u,v,y) in C2xC2xC1: f(s,t)=f(u,v), (s,t)!=(u,v) }
    # where
    #   Cn denotes n-dimensional complex space,
    #   f(s,t)=(p1:p2:p3:p4:p5), and
    #   f(u,v)=(q1:q2:q3:q4:q5).
    #
    # Notice that the variable y is used to encode the inequality.
    # We obtain the preimage of a singular component by eliminating
    # the variables {u,v,y} in getFactorsElimId(...).
    #
    return [
        p1*q2-p2*q1, p1*q3-p3*q1, p1*q4-p4*q1, p1*q5-p5*q1,
        p2*q3-p3*q2, p2*q4-p4*q2, p2*q5-p5*q2,
        p3*q4-p4*q3, p3*q5-p5*q3,
        p4*q5-p5*q4,
        (s*v-t*u)*y-1
    ]
end proc:

getFactorsElimId := proc(polyList::list)
    description "Returns factors of elements in elemination ideal of a list of polynomials.":
    #
    # Example
    #
    # > D4a:=[1,1,0,0,1,0,1,0]:
    # > lprint(getFactorsElimId(getPolyList(D4a)));
    #
    #   [
    #       [s-RootOf(3*_Z^2+4*_Z+3), 1],
    #       [4/3+RootOf(3*_Z^2+4*_Z+3)+s, 1],
    #       [t+2/5+3/5*RootOf(3*_Z^2+4*_Z+3), 1],
    #       [t-2/5-3/5*RootOf(3*_Z^2+4*_Z+3), 1],
    #       [1/2*t+t*s-1/2 *s, 1],
    #       [3*t+t*s+1+s, 1]
    #   ]
    #
    #
    local gb;

    gb:=Groebner[Basis](polyList, lexdeg([y,u,v],[s,t])):
    gb:=remove(has, gb, [u,v,y]):

    if nops(gb)<>1 then
        lprint("Warning, unexpected number of generators in elimination ideal: ", nops(gb), gb);
        if nops(gb)=0 then
            return [[0,0]];
        end if;
    end if;

    # The ouput of evala(AFactors(...)) is of the form [u,[f1,e1],...,[fn,en]]
    # and corresponds to the factorization u*f1^e1*...*fn^en
    #
    return evala(AFactors(gb[1]))[2];
end proc:

getPreSing := proc(paramType::list, paramTypeName::string:="Unknown")
    description "Computes preimage of the singular locus of the image of the map that is defined by the parametric type input.":
    #
    #
    # Example
    #
    # > getPreSing(D5c,"D5c"):
    #
    # parametric type="D5c" value=[1, 0, 0, 0, 1, 5/2, 0, 0]
    #         bideg=[0, 1], mult=1, factor=t-RootOf(5*_Z^2+8*_Z+5), sing=[1]
    #         bideg=[0, 1], mult=1, factor=8/5+RootOf(5*_Z^2+8*_Z+5)+t, sing=[1]
    #         bideg=[0, 1], mult=1, factor=t+40/33+RootOf(33*_Z^2+40*_Z+33), sing=[1]
    #         bideg=[0, 1], mult=1, factor=t-RootOf(33*_Z^2+40*_Z+33), sing=[1]
    #         bideg=[1, 0], mult=2, factor=s-4/3-5/3*RootOf(5*_Z^2+8*_Z+5), sing=[1]
    #         bideg=[1, 0], mult=2, factor=s+4/3+5/3*RootOf(5*_Z^2+8*_Z+5), sing=[1]
    #
    #
    local fctList,i,gb;
    fctList := getFactorsElimId(getPolyList(paramType)):
    printf("parametric type=%a value=%a\n", paramTypeName, paramType);
    for i from 1 to nops(fctList) do
        gb:=Groebner[Basis]([poly,diff(poly,s),diff(poly,t)], plex(s,t)):
        printf("\tbideg=%a, sing=%a, mult=%a, factor=%a\n",[degree(fctList[i][1],s),degree(fctList[i][1],t)],gb,fctList[i][2],fctList[i][1]):
    end do;
end proc:

allParamTypes := proc()
    description "Computes preimages singular locus for each predefined parametric type.":
    getPreSing(D3a,"D3a"):
    getPreSing(D3b,"D3b"):
    getPreSing(D3c,"D3c"):
    getPreSing(D3d,"D3d"):
    getPreSing(D3e,"D3e"):
    getPreSing(D4a,"D4a"):
    getPreSing(D4b,"D4b"):
    getPreSing(D4c,"D4c"):
    getPreSing(D4d,"D4d"):
    getPreSing(D5a,"D5a"):
    getPreSing(D5b,"D5b"):
    getPreSing(D5c,"D5c"):
end proc:

# do computations for all predefined parametric types.
allParamTypes();
