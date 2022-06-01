import numpy as np
import sys
import pathlib
from load_data_ import load_data
from reho_ import preproc, correlation_pearson, plot_corr
from distributions_ import run_distributions
from statistics_ import plot_feature_distribution

active_dir = pathlib.Path(sys.argv[1])
laser_folder = pathlib.Path(sys.argv[2])
vid_folder = pathlib.Path(sys.argv[3])
if not laser_folder.is_dir() or not vid_folder.is_dir():
    print('Error!: laser/video folder not found!')
    raise Exception
# imgname_fmt = "bee-001_Cycle00001_Ch2_{:06d}.ome.tif"
sleep_correlation_fname = "correlation_awake_imgs.jpg"
awake_correlation_fname = "correlation_sleep_imgs.jpg"
model_fname = "rfc.pkl"
pair_plot_fname = 'pairplot.jpg'
feature_importance_fname = "rfc_feature_importance.jpg"
roc_curve_fname = "rfc_roc_curve.jpg"
feature_distributions_fmt = 'feature_distribution_{}.jpg'


x1 = 70
x2 = 80
y1 = 70
y2 = 80
LENGTH = 100

awake_imgs, sleep_imgs = load_data(
    laser_folder=laser_folder, vid_folder=vid_folder)
print(
    f"Shape of awake and sleep datasets: awake.shape={awake_imgs.shape}, sleep.shape={sleep_imgs.shape}")

if not active_dir.is_dir():
    active_dir.mkdir()

pathlib.os.chdir(active_dir)

folder_tosave = pathlib.Path('saved')
if not folder_tosave.is_dir():
    folder_tosave.mkdir()
np.save(folder_tosave.joinpath(
    'awake_imgs.npy').as_posix(), awake_imgs)
np.save(folder_tosave.joinpath(
    'sleep_imgs.npy').as_posix(), sleep_imgs)


awake_imgs = preproc(awake_imgs)
sleep_imgs = preproc(sleep_imgs)
print(
    f"Shape of awake and sleep datasets after perproccessing: {awake_imgs.shape}, {sleep_imgs.shape}")
corr_awake_imgs = correlation_pearson(awake_imgs, neighs=5)
corr_sleep_imgs = correlation_pearson(sleep_imgs, neighs=5)
plot_corr(corr_awake_imgs, file_name=awake_correlation_fname)
plot_corr(corr_sleep_imgs, file_name=sleep_correlation_fname)
df, rfc_model = run_distributions(awake_imgs, sleep_imgs, x1=x1, x2=x2, y1=y1, y2=y2,
                                  model_fname=folder_tosave.joinpath(model_fname).as_posix(), pair_plot_fname=pair_plot_fname,
                                  feature_importance_fname=feature_importance_fname, roc_curve_fname=roc_curve_fname)
df.to_csv(folder_tosave.joinpath('features_dataframe.csv').as_posix())
plot_feature_distribution(df, save_file_fmt=feature_distributions_fmt)
