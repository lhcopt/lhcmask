import os

def make_links(links_dict, force=False):
    for kk in links_dict.keys():
        if force:
            if os.path.exists(kk):
                os.remove(kk)
        os.symlink(links_dict[kk], kk)
