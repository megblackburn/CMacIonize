def get_tb_grid(grid,subgriddim,gridsize):
    result = np.zeros((gridsize[0],gridsize[1],gridsize[2]))

    cx = int(gridsize[0]/subgriddim[0])
    cy = int(gridsize[1]/subgriddim[1])
    cz = int(gridsize[2]/subgriddim[2])

    startchunk = 0
    endchunk = cx*cy*cz
    ix = 0
    iy = 0
    iz = 0
    while endchunk <= gridsize[0]*gridsize[1]*gridsize[2]:
        chunk = np.array(grid[startchunk:endchunk])
        result[ix:ix+cx,iy:iy+cy,iz:iz+cz] = chunk.reshape(cx,cy,cz)
        startchunk += cx*cy*cz
        endchunk += cx*cy*cz
        iz += cz
        if iz == gridsize[2]:
            iz = 0
            iy += cy
            if iy == gridsize[1]:
                iy = 0
                ix += cx

    return result
