import sfsimodels as sm
import geofound as gf
import re

# TODO: Move check_required, find_and_set to sfsimodels package and raise ModelError from there


def find_and_set(objs, error):

    pattern = r"'([A-Za-z0-9_\.-]*)'"
    m = re.findall(pattern, str(error))
    parameter = m[0]
    base_type = m[1]
    for obj in objs:
        if obj.base_type == base_type:
            fn_name = "def set_%s_" % parameter
            # search current file for setters, TODO: should be able to pass modules in
            a = open(__file__)
            lines = a.readlines()
            fn = ""
            for line in lines:
                if fn_name in line:
                    fn = line.split("(")[0]
                    fn = fn.split(" ")[-1]
                    break
            local = 1
            if local:
                globals()[fn](obj)

            else:
                getattr(gf, fn)(obj)
    pass


def set_g_mod_via_unit_weight(sl):
    sl.g_mod = 13000

def run():

    sl = sm.Soil()
    sl.unit_dry_weight = 16000
    sl.e_max = 1.0
    sl.e_min = 0.5
    sl.poissons_ratio = 0.3

    fd = sm.FoundationRaft()
    fd.width = 5  # m
    fd.length = 2  # m
    fd.depth = 0

    try:
        k_rot = gf.rotational_stiffness(sl, fd)
    except gf.DesignError as e:
        print("failed, try setting with available expressions")
        find_and_set([sl, fd], e)
        k_rot = gf.rotational_stiffness(sl, fd)
    print("k_rot: ", k_rot)


if __name__ == '__main__':
    run()
