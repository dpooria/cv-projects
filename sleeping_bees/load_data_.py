import re
import cv2
import pathlib
import numpy as np
from rich.progress import track


def load_data(laser_folder='laser_data', vid_folder='vid_data',
              laser_pattern: re.Pattern = re.compile(
        r"(bee-)(\d+)\_(Cycle)(\d+)\_(Ch)(\d+)\_(\d+)\.(ome)\.(tif)"),
        vid_pattern: re.Pattern = re.compile(
            r"(\d+)\_(sleep|awake)\_(\d+)\.mat")
):
    laser_folder = pathlib.Path(laser_folder)
    vid_folder = pathlib.Path(vid_folder)
    if not laser_folder.is_dir() or not vid_folder.is_dir():
        print('Error! no such directories for video and laser')
        raise Exception()
    laser_data = [l.relative_to(laser_folder).as_posix() for l in laser_folder.iterdir(
    ) if laser_pattern.search(l.relative_to(laser_folder).as_posix())]
    vid_data = [v.relative_to(vid_folder).as_posix() for v in vid_folder.iterdir(
    ) if vid_pattern.search(v.relative_to(vid_folder).as_posix())]
    laser_cycles = [int(laser_pattern.sub(r'\4', l)) for l in laser_data]
    vid_cycle_label = {int(vid_pattern.sub(
        r'\1', v)): vid_pattern.sub(r'\2', v) for v in vid_data}
    awake_imgs = []
    sleep_imgs = []
    laser_folder_path = pathlib.Path(laser_folder)
    for im_name, cycle in track(zip(laser_data, laser_cycles),
                                f'Loading laser images from the folder {laser_folder.as_posix()}',
                                total=len(laser_data)):
        tmp = cv2.imread(str(laser_folder_path.joinpath(im_name)),
                         cv2.IMREAD_GRAYSCALE)
        if vid_cycle_label[cycle] == 'sleep':
            sleep_imgs.append(tmp)
        elif vid_cycle_label[cycle] == 'awake':
            awake_imgs.append(tmp)
    awake_imgs = np.asarray(awake_imgs)
    sleep_imgs = np.asarray(sleep_imgs)

    
    return awake_imgs, sleep_imgs
