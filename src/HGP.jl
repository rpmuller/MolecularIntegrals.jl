export coulomb_hgp

function coulomb_hgp(a::PGBF,b::PGBF,c::PGBF,d::PGBF)
    return a.norm*b.norm*c.norm*d.norm*hrr(a.expn,a.x,a.y,a.z,a.I,a.J,a.K,
        b.expn,b.x,b.y,b.z,b.I,b.J,b.K,
        c.expn,c.x,c.y,c.z,c.I,c.J,c.K,
        d.expn,d.x,d.y,d.z,d.I,d.J,d.K)
end
coulomb_hgp(a::CGBF,b::CGBF,c::CGBF,d::CGBF) = contract(coulomb_hgp,a,b,c,d)

function hrr(aexpn,ax,ay,az,aI,aJ,aK,
    bexpn,bx,by,bz,bI,bJ,bK,
    cexpn,cx,cy,cz,cI,cJ,cK,
    dexpn,dx,dy,dz,dI,dJ,dK)
    if bI > 0
        return hrr(aexpn,ax,ay,az,aI+1,aJ,aK,bexpn,bx,by,bz,bI-1,bJ,bK,
            cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK) +
        (ax-bx)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI-1,bJ,bK,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK)
    elseif bJ > 0
        return hrr(aexpn,ax,ay,az,aI,aJ+1,aK,bexpn,bx,by,bz,bI,bJ-1,bK,
            cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK) +
            (ay-by)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ-1,bK,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK)
    elseif bK > 0
        return hrr(aexpn,ax,ay,az,aI,aJ,aK+1,bexpn,bx,by,bz,bI,bJ,bK-1,
            cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK) +
            (az-bz)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK-1,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK)
    elseif dI > 0
        return hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
            cexpn,cx,cy,cz,cI+1,cJ,cK,dexpn,dx,dy,dz,dI-1,dJ,dK) +
            (cx-dx)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI-1,dJ,dK)
    elseif dJ > 0
        return hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
            cexpn,cx,cy,cz,cI,cJ+1,cK,dexpn,dx,dy,dz,dI,dJ-1,dK) +
            (cy-dy)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ-1,dK)
    elseif dK > 0
        return hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
            cexpn,cx,cy,cz,cI,cJ,cK+1,dexpn,dx,dy,dz,dI,dJ,dK-1) +
            (cz-dz)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK-1)
    end
    return vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
               cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,0)
end


function vrr(aexpn,ax,ay,az,aI,aJ,aK,
        bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI,cJ,cK,
        dexpn,dx,dy,dz,m)
    px,py,pz = gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
    qx,qy,qz = gaussian_product_center(cexpn,cx,cy,cz,dexpn,dx,dy,dz)
    zeta,eta = aexpn+bexpn,cexpn+dexpn
    wx,wy,wz = gaussian_product_center(zeta,px,py,pz,eta,qx,qy,qz)
    
    val = 0
    if cK>0
        val = (qz-cz)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
            cexpn,cx,cy,cz,cI,cJ,cK-1,dexpn,dx,dy,dz,m) +
            (wz-qz)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK-1,dexpn,dx,dy,dz,m+1)
        if cK>1
            val += 0.5*(cK-1)/eta*(
                vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                    cexpn,cx,cy,cz,cI,cJ,cK-2,dexpn,dx,dy,dz,m) -
            zeta/(zeta+eta)*
                vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                    cexpn,cx,cy,cz,cI,cJ,cK-2,dexpn,dx,dy,dz,m+1) )
        end
        if aK>0
            val += 0.5*aK/(zeta+eta)*
                vrr(aexpn,ax,ay,az,aI,aJ,aK-1,bexpn,bx,by,bz,
                    cexpn,cx,cy,cz,cI,cJ,cK-1,dexpn,dx,dy,dz,m+1)
        end
        return val
    elseif cJ>0
        val = (qy-cy)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
            cexpn,cx,cy,cz,cI,cJ-1,cK,dexpn,dx,dy,dz,m) +
        (wy-qy)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ-1,cK-1,dexpn,dx,dy,dz,m+1)
        #println("val4=$val")
        if cJ>1
            val += 0.5*(cJ-1)/eta*(
            vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ-2,cK,dexpn,dx,dy,dz,m) -
            zeta/(zeta+eta)*
            vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ-2,cK,dexpn,dx,dy,dz,m+1)
            )
        #println("val5=$val")
        end
        if aJ>0
            val += 0.5*aJ/(zeta+eta)*
            vrr(aexpn,ax,ay,az,aI,aJ-1,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ-1,cK,dexpn,dx,dy,dz,m+1)
        end
        #println("val6=$val")
        return val
    elseif cI>0
        val = (qx-cx)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
            cexpn,cx,cy,cz,cI-1,cJ,cK,dexpn,dx,dy,dz,m) +
        (wx-qx)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI-1,cJ,cK-1,dexpn,dx,dy,dz,m+1)
        #println("val7=$val")
        if cI>1
            val += 0.5*(cI-1)/eta*(
            vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI-2,cJ,cK,dexpn,dx,dy,dz,m) -
            zeta/(zeta+eta)*
            vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI-2,cJ,cK,dexpn,dx,dy,dz,m+1)
            )
        end
        #println("val8=$val")
        if aI>0
            val += 0.5*aI/(zeta+eta)*
            vrr(aexpn,ax,ay,az,aI-1,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI-1,cJ,cK,dexpn,dx,dy,dz,m+1)
        end
        #println("val9=$val")
        return val
    elseif aK>0
        val = (pz-az)*vrr(aexpn,ax,ay,az,aI,aJ,aK-1,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m) +
        (wz-pz)*vrr(aexpn,ax,ay,az,aI,aJ,aK-1,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m+1)
        #println("val10=$val")
        if aK>1
            val += 0.5*(aK-1)/zeta*(
            vrr(aexpn,ax,ay,az,aI,aJ,aK-2,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI-1,cJ,cK,dexpn,dx,dy,dz,m) -
            eta/(zeta+eta)*
            vrr(aexpn,ax,ay,az,aI,aJ,aK-2,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI-1,cJ,cK,dexpn,dx,dy,dz,m+1)
            )
        end
        #println("val11=$val")
        return val
    elseif aJ>0
        val = (py-ay)*vrr(aexpn,ax,ay,az,aI,aJ-1,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m)+
        (wy-py)*vrr(aexpn,ax,ay,az,aI,aJ-1,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m+1)
        #println("val12=$val")
        if aJ>1
            val += 0.5*(aJ-1)/zeta*(
            vrr(aexpn,ax,ay,az,aI,aJ-2,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m) -
            eta/(zeta+eta)*
            vrr(aexpn,ax,ay,az,aI,aJ-2,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m+1)
            )
        end
        #println("val13=$val")
        return val
    elseif aI>0
        val = (px-ax)*vrr(aexpn,ax,ay,az,aI-1,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m) +
        (wx-px)*vrr(aexpn,ax,ay,az,aI-1,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m+1)
        #println("val14=$val")
        if aI>1
            val += 0.5*(aI-1)/zeta*(
            vrr(aexpn,ax,ay,az,aI-2,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m) -
            eta/(zeta+eta)*
            vrr(aexpn,ax,ay,az,aI-2,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m+1)
            )
        end
        #println("val15=$val")
        return val
    end

    rab2 = dist2(ax-bx,ay-by,az-bz)
    rcd2 = dist2(cx-dx,cy-dy,cz-dz)
    rpq2 = dist2(px-qx,py-qy,pz-qz)
    T = zeta*eta/(zeta+eta)*rpq2
    Kab = sqrt(2)pi^1.25/zeta*exp(-aexpn*bexpn*rab2/zeta)
    Kcd = sqrt(2)pi^1.25/eta*exp(-cexpn*dexpn*rcd2/eta)
    #println("rab2=$rab2,rcd2=$rcd2,rpq2=$rpq2,T=$T,Kab=$Kab,Kcd=$Kcd")
    return Kab*Kcd/sqrt(zeta+eta)*Fgamma(m,T)
end

