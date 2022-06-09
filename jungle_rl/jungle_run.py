
from itertools import cycle
import tensorflow as tf
import numpy as np
from tf_agents.environments import TFPyEnvironment

import time
import pathlib
import curses
from functools import partial

from agent import IMAgent
from jenv import JungleEnv


grid = (20, 20)
L = 8
N = 4


def rescale(x, y, w, h):
    if x == 0:
        rx = 0
    elif x == grid[0] - 1:
        rx = w - 2
    else:
        rx = int((x + 1) / grid[0] * w) - 2

    if y == 0:
        ry = 0
    elif y == grid[1] - 1:
        ry = h - 2
    else:
        ry = int((y + 1) / grid[1] * h) - 2

    return [ry, rx]


def jungle_action_fn(agent, action):
    return {'action': action, 'agent': agent}


policy_dir = pathlib.Path('saved_policies-1').joinpath('saved_policies')
# policy_path = [policy_dir.joinpath(f'agent_{i+2}') for i in range(N)]
policy_path = [policy_dir.joinpath(
    f'agent_2')] * 2 + [policy_dir.joinpath(f'agent_4')] * 2

policy_list = [tf.saved_model.load(str(pp)) for pp in policy_path]
eval_jungle_env = JungleEnv(grid=grid, N=N, L=L)
tfeval_jungle_env = TFPyEnvironment(eval_jungle_env)
agent_list = [IMAgent(tfeval_jungle_env, action_spec=tfeval_jungle_env.action_spec()['action'],
                      action_fn=partial(jungle_action_fn, i + 2),
                      name=f'agent{i+2}', id=i+2) for i in range(N)]


# initial screen
sc = curses.initscr()
h, w = sc.getmaxyx()
win = curses.newwin(h, w, 0, 0)


win.keypad(0)
curses.curs_set(0)

#define colors (256 based)
RED = 9
GREEN = 10
YELLOW = 11
BLUE = 19
WHITE = 15


agent_colors = [GREEN, YELLOW, BLUE, WHITE]
food_color = RED
colors = curses.can_change_color()
if colors:
    curses.start_color()
    curses.use_default_colors()
    curses.init_pair(RED, curses.COLOR_RED, -1)
    curses.init_pair(GREEN, curses.COLOR_GREEN, -1)
    curses.init_pair(YELLOW, curses.COLOR_YELLOW, -1)
    curses.init_pair(BLUE, curses.COLOR_BLUE, -1)
    curses.init_pair(WHITE, curses.COLOR_WHITE, -1)


def next_step(agent, policy):
    agent.act_eval(policy)
    ts = tfeval_jungle_env._current_time_step()
    observation = ts.observation.numpy()[0, :, :]
    state = observation[:, :-1]
    rewards = observation[:N, -1]
    return state, ts.is_last(), rewards


def run():
    ts = tfeval_jungle_env.reset()
    for a in agent_list:
        a.reset()

    random_priority = np.random.permutation(range(N))
    agent_rand = list(np.array(agent_list)[random_priority])
    policy_rand = list(np.array(policy_list)[random_priority])
    agent_cycle = cycle(zip(agent_rand, policy_rand))

    agent, policy = next(agent_cycle)
    state, is_last, rewards = next_step(agent, policy)

    # foods
    agent_positions = [np.where(state == i+2) for i in range(N)]
    fl = np.where(state == 1)
    food_positions = [(x, y) for x, y in zip(fl[0], fl[1])]

    iter = 0

    ts = tfeval_jungle_env._current_time_step()
    rewards = ts.observation.numpy()[0, :N, -1]
    win.addstr(
        0, 0, f'step: {iter:3d} agent2_reward: {rewards[0]:3d}, agent3_reward: {rewards[1]:3d}, agent4_reward: {rewards[2]:3d}, agent5_reward: {rewards[3]:3d}')

    while not is_last:
        win.refresh()
        time.sleep(0.1)
        prev_agent_positions = np.copy(agent_positions)
        prev_food_positions = np.copy(food_positions)
        prev_rewards = np.copy(rewards)

        iter += 1
        win.border(0)
        win.timeout(100)
        agent, policy = next(agent_cycle)
        state, is_last, rewards = next_step(agent, policy)
        win.addstr(
            0, 0, f'step: {iter:3d} agent2_reward: {rewards[0]:3d}, agent3_reward: {rewards[1]:3d}, agent4_reward: {rewards[2]:3d}, agent5_reward: {rewards[3]:3d}')

        # foods
        agent_positions = [np.where(state == i+2) for i in range(N)]
        fl = np.where(state == 1)
        food_positions = [(x, y) for x, y in zip(fl[0], fl[1])]
        for f, fp in zip(food_positions, prev_food_positions):
            rfp = rescale(fp[0], fp[1], w, h)
            rf = rescale(f[0], f[1], w, h)
            win.addch(rfp[0], rfp[1], ' ')
            if colors:
                win.addch(rf[0], rf[1], curses.ACS_DIAMOND,
                          curses.color_pair(food_color))
            else:
                win.addch(rf[0], rf[1], curses.ACS_DIAMOND)

        for i, (a, ap) in enumerate(zip(agent_positions, prev_agent_positions)):
            rap = rescale(ap[0][0], ap[1][0], w, h)
            win.addch(rap[0], rap[1], ' ')

            ra = rescale(a[0][0], a[1][0], w, h)
            if colors:
                win.addch(ra[0], ra[1], curses.ACS_BLOCK,
                          curses.color_pair(agent_colors[i]))
            else:
                win.addch(ra[0], ra[1], str(i))

    win.clear()
    sc.addstr(
        10, 0, f'finished!: agent2_reward: {prev_rewards[0]:3d}, agent3_reward: {prev_rewards[1]:3d},'
        f' agent4_reward: {prev_rewards[2]:3d}, agent5_reward: {prev_rewards[3]:3d}')

    sc.refresh()
    time.sleep(5)


try:
    while True:
        win.addstr(h//2 - 1, w//2 - 10, 'press ENTER or ESC to exit')
        win.addstr(h//2, w//2 - 10, 'press space to run!')
        next_key = win.getch()

        if next_key == 10 or next_key == 27:
            curses.endwin()
            break
        elif next_key != -1:
            win.clear()
            run()

except:
    curses.endwin()
