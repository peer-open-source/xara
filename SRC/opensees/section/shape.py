import numpy as np



def create_rectpatch(d, b, yoff, zoff):
    return np.array([
            [yoff - d / 2, zoff - b / 2],
            [yoff + d / 2, zoff - b / 2],
            [yoff + d / 2, zoff + b / 2],
            [yoff - d / 2, zoff + b / 2],
            # [yoff - d / 2, zoff - b / 2],
    ]).T


def rectangle2fiber(yz, int_typ, ny, nz):
    y = np.linspace(yz[0, 0], yz[1, 0], ny + 1)
    z = np.linspace(yz[0, 1], yz[1, 1], nz + 1)
    dy = np.diff(y)
    dz = np.diff(z)
    y = y[:-1] + dy / 2
    z = z[:-1] + dz / 2
    y, z = np.meshgrid(y, z)
    yfib = y.ravel()
    zfib = z.ravel()
    wfib = dy * dz
    return yfib, zfib, wfib


def create_ipmesh4wfshape(sec_data, ndm=3):
    d  = sec_data['d']
    bf = sec_data['bf']
    tf = sec_data['tf']
    tw = sec_data['tw']

    nft = sec_data['nft']
    nwl = sec_data['nwl']
    nfl = sec_data.get('nfl', 1, ) #sec_data['nft'])
    nwt = sec_data.get('nwt', 1, ) #sec_data['nwl'])

    int_typ = sec_data.get("IntTyp", None)
    flg_opt = sec_data.get('FlgOpt', True)
    web_opt = sec_data.get('WebOpt', True)

    yoff = ( d - tf) / 2
    zoff = (bf + tw) / 4

    dw = d - 2 * tf
    bi = bf - tw
    if flg_opt:
        yz = np.vstack([create_rectpatch(tf, bf,  yoff, 0),
                        create_rectpatch(tf, bf, -yoff, 0)])
        outline = [yz.T]
        yz = create_rectpatch(dw, tw, 0, 0)
        outline.append(yz.T)
    else:
        if web_opt:
            yz = np.vstack([create_rectpatch(tf, bi / 2,  yoff,  zoff),
                            create_rectpatch(tf, bi / 2,  yoff, -zoff),
                            create_rectpatch(tf, bi / 2, -yoff,  zoff),
                            create_rectpatch(tf, bi / 2, -yoff, -zoff)])
            outline = [yz.T]
            yz = create_rectpatch(d, tw, 0, 0)
            outline.append(yz.T)
        else:
            yz = np.vstack([create_rectpatch(tf, bi / 2, yoff, zoff),
                            create_rectpatch(tf, bi / 2, yoff, -zoff),
                            create_rectpatch(tf, bi / 2, -yoff, zoff),
                            create_rectpatch(tf, bi / 2, -yoff, -zoff)])
            outline = [yz.T]
            yz = create_rectpatch(dw, tw, 0, 0)
            outline.append(yz.T)
            dip = tf
            bip = tw
            yz = np.vstack([create_rectpatch(dip, bip, yoff, 0),
                            create_rectpatch(dip, bip, -yoff, 0)])
            outline.append(yz.T)

    mesh = [{'ny': [nft, nft, nft, nft], 
             'nz': [nfl, nfl, nfl, nfl]}, 
            {'ny': [nwl, nwl], 
             'nz': [nwt, nwt]}]

    if not flg_opt and not web_opt:
        niy = sec_data['niy']
        niz = sec_data['niz']
        mesh.append({'ny': [niy, niy, niy, niy], 
                     'nz': [niz, niz, niz, niz]})

    no = len(outline)
    fiby = [None] * no
    fibz = [None] * no
    fiba = [None] * no
    for oc in range(no):
        ny   = mesh[oc]["ny"][0]
        nz   = mesh[oc]["nz"][0]
        patch_cor = outline[oc]
        nip = patch_cor.shape[1]
        yfib = np.zeros((ny * nz, nip))
        zfib = np.zeros((ny * nz, nip))
        wfib = np.zeros((ny * nz, nip))
        for ip in range(nip):
            ny = mesh[oc]['ny'][ip]
            nz = mesh[oc]['nz'][ip]
            yz = np.reshape(patch_cor[:, ip], (2, 2))
            yfib[:, ip], zfib[:, ip], wfib[:, ip] = rectangle2fiber(yz.T, int_typ, ny, nz)
        fiby[oc] = yfib.ravel()
        fibz[oc] = zfib.ravel()
        fiba[oc] = wfib.ravel()

    yfib   = np.concatenate(fiby)
    zfib   = np.concatenate(fibz)
    wfib   = np.concatenate(fiba)
    mat_id = np.ones(len(yfib))

    return -zfib, yfib,  wfib, mat_id


