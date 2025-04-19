
class Analysis:
    def __init__(self, model):
        self._model = model
        self._test = None

        self._integrator = None

    def integrator(self, *args):
        self._integrator = args

    def increment(self, factor=1.0):
        self._model.integrator("LoadControl", factor)
        self._model.analysis("Static")
        self._model.analyze(1, operation='increment')

    def iteration(self):
#       self._model.test("FixedNumIter", 1, 1)
        niter = 1
        tol = 1e-12
        self._model.test("EnergyIncr", tol, niter, 5)
        self._model.analysis("Static")
        self._model.analyze(1, operation='iteration')

def _buckle_factor(boundary, phi=0):
    import math
    if boundary == "pin-pin":
        return math.pi

    if boundary == "fix-slide":
        return math.pi

    if boundary == "fix-fix":
        return 2*math.pi

    if boundary == "fix-pin":
        import scipy.optimize
        f = lambda x: math.tan(x) - x/(1 + x**2*phi/12)
        sol = scipy.optimize.root_scalar(f, x0=0.7, bracket=(math.pi, 1.45*math.pi))
        if sol.converged:
            return sol.root

    if boundary == "fix-free":
        return math.pi/2

    if boundary == "pin-slide":
        return math.pi/2

def _fix_node(model, node, type):
    ndf = model.getNDF()
    reactions = [0 for _ in range(ndf)]

    long, tran, vert = range(3)
    if ndf == 6:
        bend = vert+3
        # always fix out-of-plane rotation, which spins
        # about the transverse DOF
        reactions[tran+3] = 1
        reactions[long+3] = 1 # torsion
        reactions[2] = 1
    else:
        bend = vert


    if node > 1:
        vert = 0
    else:
        vert = 1

    if type == "fix":
        reactions[tran] = 1
        reactions[long] = 0 if node > 1 else 1
        reactions[bend] = 1

    elif type == "pin":
        reactions[tran] = 1
        reactions[long] = 0 if node > 1 else 1
        reactions[bend] = 0

    elif type == "slide":
        reactions[tran] = 0
        reactions[long] = 0 if node > 1 else 1
        reactions[bend] = 1

    elif type == "free":
        pass

    model.fix(node, *reactions)
    return reactions


def create_prism(length:    float,
                 element:   str,
                 section:   dict,
                 boundary:  tuple,
                 geometry:  str = None,
                 transform: str = None,
                 divisions: int = 1,
                 rotation = None, 
                 ndm=3):
    import opensees.openseespy as ops

    L  = length

    # Number of elements discretizing the column
    ne = divisions

    elem_type  = element
    geom_type  = transform

    # Number of integration points along each element
    nIP = 5
    nn = ne + 1

    model = ops.Model(ndm=ndm)

    for i in range(1, nn+1):
        x = (i-1)/float(ne)*L
        if ndm == 3:
            location = (x, 0.0, 0.0)
        else:
            location = (x, 0.0)

        if rotation is not None:
            location = tuple(rotation@location)

        model.node(i, location)

        model.mass(i, *[1.0]*model.getNDF())


    # Define boundary conditions
    if isinstance(boundary[0], str):
        _fix_node(model,  1, boundary[0])
    else:
        model.fix(1, boundary[0])

    if isinstance(boundary[1], str):
        _fix_node(model, nn, boundary[1])
    else:
        model.fix(nn, boundary[1])

    #
    # Define cross-section 
    #
    sec_tag = 1
    properties = []
    for k,v in section.items():
        properties.append("-" + k)
        properties.append(v)

    model.section("FrameElastic", sec_tag, *properties)

    # Define geometric transformation
    geo_tag = 1
    if ndm == 3:
        vector = (0,  0, 1)
#           vector = (0,  -1, 0)
        if rotation is not None:
            vector = tuple(map(float, rotation@vector))
        # for |_ 
    else:
        vector = ()

    model.geomTransf(geom_type, geo_tag, *vector)

    # Define elements
    for i in range(1, ne+1):
        if geometry is None or geometry == "Linear" or "Exact" in elem_type:
            model.element(elem_type, i, (i, i+1),
                        section=sec_tag,
                        transform=geo_tag)
        else:
            model.element(elem_type, i, (i, i+1),
                        section=sec_tag,
                        order={"Linear": 0, "delta": 1}[geometry],
                        transform=geo_tag)

    return model

class Prism:
    def __init__(self,
                 length:    float,
                 element:   str,
                 section:   dict,
                 boundary:  tuple,
                 geometry:  str = None,
                 transform: str = None,
                 divisions: int = 1,
                 rotation = None,
                 shear = True):

        self.length    = length
        self.element   = element
        self.section   = section
        self.boundary  = boundary
        self.geometry  = geometry
        self.transform = transform
        self.divisions = divisions
        self.use_shear = shear
        self.rotation = rotation

    def create_model(self, ndm=3, file=None):
        import opensees.openseespy as ops

        L  = self.length
        boundary = self.boundary

        # Number of elements discretizing the column
        ne = self.divisions

        elem_type  = self.element
        geom_type  = self.transform

        # Number of integration points along each element
        nIP = 5
        nn = ne + 1
        if file is not None:
            file = open(file, "w+")
        model = ops.Model(ndm=ndm, echo_file=file)

        for i in range(1, nn+1):
            x = (i-1)/float(ne)*L
            if ndm == 3:
                location = (x, 0.0, 0.0)
            else:
                location = (x, 0.0)

            if self.rotation is not None:
                location = tuple(self.rotation@location)

            model.node(i, location)

            model.mass(i, *[1.0]*model.getNDF())


        # Define boundary conditions
        if isinstance(boundary[0], str):
            _fix_node(model,  1, boundary[0])
        else:
            model.fix(1, boundary[0])

        if isinstance(boundary[1], str):
            _fix_node(model, nn, boundary[1])
        else:
            model.fix(nn, boundary[1])

        #
        # Define cross-section 
        #
        sec_tag = 1
        properties = []
        for k,v in self.section.items():
            properties.append("-" + k)
            properties.append(v)

        model.section("ElasticFrame", sec_tag, *properties)

        # Define geometric transformation
        geo_tag = 1
        if ndm == 3:
            vector = (0,  0, 1)
#           vector = (0,  -1, 0)
            if self.rotation is not None:
                vector = tuple(map(float, self.rotation@vector))
            # for |_ 
        else:
            vector = ()

        model.geomTransf(geom_type, geo_tag, *vector)

        # Define elements
        for i in range(1, ne+1):
            if self.geometry is None or self.geometry == "Linear" or "Exact" in elem_type:
                model.element(elem_type, i, (i, i+1),
                            section=sec_tag,
                            transform=geo_tag,
#                           iter=(2, 1e-12),
                            shear=int(self.use_shear))
            else:
                model.element(elem_type, i, (i, i+1),
                            section=sec_tag,
                            order={"Linear": 0, "delta": 1}[self.geometry],
                            transform=geo_tag,
                            shear=int(self.use_shear))

        return model

    def buckle(self, model):
        # Define loads
        if self.use_shear:
            phi = 12*E*I/(Ay*G*L**2)
            lam = _buckle_factor(boundary, phi)
            kL = L/lam
            euler_load = E*I/kL**2  / (1 + lam**2*phi/12)
        else:
            kL = L/_buckle_factor(boundary)
            euler_load = E*I/kL**2

        model.pattern("Plain", 1, "Linear")
        if self.ndm == 2:
            load = (-euler_load, 0.0, 0.0)
        else:
            load = (-euler_load, 0.0, 0.0, 0, 0, 0)

        model.load(nn, *load, pattern=1)

        return model


def node_average(model, response, ndm=2, keys=None):
    """
    Return a dictionary of stress components per node.
    For 2D: sxx, syy, sxy
    For 3D: sxx, syy, szz, sxy, syz, sxz
    Summed values from each element's node contribution are averaged at the end.

    Args:
        model: an OpenSees model-like object with the following interface:
               - getEleTags()
               - getEleClassTags(eleTag)
               - getNodeTags()
               - eleNodes(eleTag)
               - eleResponse(eleTag, str)
        response (str): a label if needed for e.g. 'stress', 'stressAtNodes'
        ndm (int): number of dimensions (2 or 3)

    Returns:
        dict[nodeTag] = {
            "xx": float, 
            "yy": float,
            (optionally in 2D) "xy": float,
            (optionally in 3D) "zz", "yz", "xz",
            ... possibly derived quantities ...
        }
    """

    # Helper: pick which elements are 2D vs 3D, if needed
    # For example, you might have separate function(s) or logic:
    #    ele_tags = _stress_2d_ele_tags_only(model, model.getEleTags())
    # or for 3D:
    #    ele_tags = _stress_3d_ele_tags_only(model, model.getEleTags())
    #
    # Below, we simply re-use the user's helper for 2D, but you'll have to define
    # your own for 3D or unify them if your model can contain both:
    ele_tags_all = model.getEleTags()
    if ndm == 2:
        ele_tags = _stress_2d_ele_tags_only(model, ele_tags_all)
    else:
        # Placeholder; you must implement a 3D filter or skip if you already know
        # the elements are 3D:
        ele_tags = ele_tags_all

    # For a typical 2D element with 3 stress components: sxx, syy, sxy
    # For a typical 3D element with 6 stress components: sxx, syy, szz, sxy, syz, sxz
    # You may need to adapt this depending on your element type(s).
    if ndm == 2:
        keys = "sxx", "syy", "sxy"
    else:
        keys = "sxx", "syy", "szz", "sxy", "syz", "sxz"


    node_tags = model.getNodeTags()

    # Initialize a dictionary to hold stress sums and a count of how many
    # times each node is encountered (to do averaging).
    sig_dict = {node: {"_count": 0, **{key: 0.0 for key in keys}} for node in node_tags}

    nrc = len(keys)

    # We assume each element’s node-level stress output is available:
    # e.g. model.eleResponse(elem_tag, 'stressAtNodes')
    # and yields a flat array (with length = n_nodes_in_element * nrc).
    # Then we reshape it to (nen, nrc).
    for elem in ele_tags:
        ele_node_tags = model.eleNodes(elem)
        nen = len(ele_node_tags)

        sig_per_node = np.reshape(model.eleResponse(elem, response), (nen, nrc))

        # Accumulate into sig_dict
        for i, node in enumerate(ele_node_tags):
            for j, key in enumerate(keys):
                sig_dict[node][key] += sig_per_node[i][j]

            sig_dict[node]["_count"] += 1

    # Now do averaging (if a node belongs to multiple elements).
    for node in node_tags:
        c = sig_dict[node]["_count"] or 1.0

        for key in keys:
            sig_dict[node][key] /= c


    # Clean up
    for node in node_tags:
        sig_dict[node].pop("_count", None)

    return sig_dict


def find_node(model, x=None, y=None, z=None, tol=1e-8):
    """
    Return the tag of the *first* node that matches the specified coordinates
    within a given tolerance. If no match is found, returns None.

    Parameters
    ----------
    model : opensees.openseespy.Model
        The model object containing nodes.
    x : float, optional
        X-coordinate to match (if None, no X match is required).
    y : float, optional
        Y-coordinate to match (if None, no Y match is required).
    tol : float
        Matching tolerance for coordinates.

    Returns
    -------
    int or None
        The tag of the first matching node, or None if not found.
    """
    # Get the list of all node tags in the model
    node_tags = model.getNodeTags()

    for tag in node_tags:
        coords = model.nodeCoord(tag)  # [xCoord, yCoord, ...]
        # Check requested coordinates
        if x is not None and abs(coords[0] - x) > tol:
            continue

        if y is not None and abs(coords[1] - y) > tol:
            continue

        if z is not None and abs(coords[2] - z) > tol:
            continue

        # If we've passed all checks, this node matches
        return tag

    # If no node passes the test, return None
    return None


def find_nodes(model, x=None, y=None, z=None, tol=1e-8):
    """
    Return an iterator of node tags that match the specified coordinates
    within a given tolerance.

    Parameters
    ----------
    model : opensees.openseespy.Model
        The model object containing nodes.
    x : float, optional
        X-coordinate to match (if None, no X match is required).
    y : float, optional
        Y-coordinate to match (if None, no Y match is required).
    tol : float
        Matching tolerance for coordinates.

    Yields
    ------
    int
        The tag of each matching node.
    """
    node_tags = model.getNodeTags()
    for tag in node_tags:
        coords = model.nodeCoord(tag)  # [xCoord, yCoord, ...]
        # Check X coordinate if requested
        if x is not None and abs(coords[0] - x) > tol:
            continue
        # Check Y coordinate if requested
        if y is not None and abs(coords[1] - y) > tol:
            continue
        # Check Z coordinate if requested
        if z is not None and abs(coords[2] - z) > tol:
            continue
        # If we've passed all checks, this node matches
        yield tag

