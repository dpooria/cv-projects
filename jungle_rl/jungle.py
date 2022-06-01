
from agent import IMAgent
from jenv import JungleEnv
from tf_agents.utils import common
from tf_agents.environments import TFPyEnvironment
from functools import partial
from itertools import cycle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()


def print_jungle(state):
    print(state)


# test environment

# jungle_env = JungleEnv()
# ts = jungle_env.reset()

# print_jungle(ts.observation)

# agents = cycle(range(2, 6))
# print(ts.is_last())
# while not np.any(ts.is_last()):
#     agent = next(agents)
#     action = {
#         'action': np.asarray(np.random.randint(0, 8)),
#         'agent': agent
#     }
#     ts = jungle_env.step(action)
#     print('Player:', agent, 'Action:', action['action'],
#           'Reward:', ts.reward, 'Board:')
#     print_jungle(ts.observation)

jungle_env = JungleEnv()
tf_jungle_env = TFPyEnvironment(jungle_env)
# # tf_jungle_env = jungle_env


def jungle_action_fn(agent, action):
    return {'action': action, 'agent': agent}


# agent_id2 = IMAgent(
#     tf_jungle_env,
#     action_spec=tf_jungle_env.action_spec()['action'],
#     action_fn=partial(jungle_action_fn, 2),
#     name='agent2', id=2
# )

# agent_id3 = IMAgent(
#     tf_jungle_env,
#     action_spec=tf_jungle_env.action_spec()['action'],
#     action_fn=partial(jungle_action_fn, 3),
#     name='agent3', id=3
# )

# agent_id4 = IMAgent(
#     tf_jungle_env,
#     action_spec=tf_jungle_env.action_spec()['action'],
#     action_fn=partial(jungle_action_fn, 4),
#     name='agent4', id=4
# )

# agent_id5 = IMAgent(
#     tf_jungle_env,
#     action_spec=tf_jungle_env.action_spec()['action'],
#     action_fn=partial(jungle_action_fn, 5),
#     name='agent5', id=5
# )

# agent_list = [agent_id2, agent_id3, agent_id4, agent_id5]

# ts = tf_jungle_env.reset()

# agents_qeue = cycle(np.random.permutation(range(2, 6)))

# ts = tf_jungle_env.current_time_step()
# state = ts.observation
# print('Random start board:')
# print_jungle(state.numpy())


# counter = 0
# while not np.any(ts.is_last()):
#     counter += 1
#     agent = agent_list[next(agents_qeue) - 2]
#     agent.act()
#     ts = tf_jungle_env.current_time_step()
#     print(f"iteration #{counter}")
#     print(f"ts_reward: {ts.reward}")
#     print(f'Agent: {agent.name}, Reward: {agent.reward}')
#     print("ts.observation:")
#     print_jungle(ts.observation.numpy())
#     input('next?')


num_iterations = 2000
initial_collect_episodes = 100
episodes_per_iteration = 10
train_steps_per_iteration = 1
training_batch_size = 512
training_num_steps = 2
replay_buffer_size = 3 * episodes_per_iteration * 9
learning_rate = 1e-2
plot_interval = 50
iteration = 1
games = []
loss_infos = []

agent_id2 = IMAgent(
    tf_jungle_env,
    action_spec=tf_jungle_env.action_spec()['action'],
    action_fn=partial(jungle_action_fn, 2),
    name='agent2', id=2, learning_rate=learning_rate,
    training_batch_size=training_batch_size,
    training_num_steps=training_num_steps,
    replay_buffer_max_length=replay_buffer_size,
    td_errors_loss_fn=common.element_wise_squared_loss
)


agent_id3 = IMAgent(
    tf_jungle_env,
    action_spec=tf_jungle_env.action_spec()['action'],
    action_fn=partial(jungle_action_fn, 3),
    name='agent3', id=3, learning_rate=learning_rate,
    training_batch_size=training_batch_size,
    training_num_steps=training_num_steps,
    replay_buffer_max_length=replay_buffer_size,
    td_errors_loss_fn=common.element_wise_squared_loss
)


agent_id4 = IMAgent(
    tf_jungle_env,
    action_spec=tf_jungle_env.action_spec()['action'],
    action_fn=partial(jungle_action_fn, 4),
    name='agent4', id=4, learning_rate=learning_rate,
    training_batch_size=training_batch_size,
    training_num_steps=training_num_steps,
    replay_buffer_max_length=replay_buffer_size,
    td_errors_loss_fn=common.element_wise_squared_loss
)


agent_id5 = IMAgent(
    tf_jungle_env,
    action_spec=tf_jungle_env.action_spec()['action'],
    action_fn=partial(jungle_action_fn, 5),
    name='agent5', id=5, learning_rate=learning_rate,
    training_batch_size=training_batch_size,
    training_num_steps=training_num_steps,
    replay_buffer_max_length=replay_buffer_size,
    td_errors_loss_fn=common.element_wise_squared_loss
)


agent_list = [agent_id2, agent_id3, agent_id4, agent_id5]


def training_episode(tf_jungle_env, agent_list):
    ts = tf_jungle_env.reset()
    agent_cycle = cycle(
        np.array(agent_list)[np.random.permutation(range(len(agent_list)))])
    for a in agent_list:
        a.reset()
    time_steps = []
    while not ts.is_last():
        player = next(agent_cycle)
        player.act(collect=True)
        ts = tf_jungle_env.current_time_step()
        time_steps.append(ts)
    return time_steps


print('Collecting Initial Training Sample...')
for _ in range(initial_collect_episodes):
    training_episode(tf_jungle_env, agent_list)
print('Samples collected')


def collect_training_data():
    for game in range(episodes_per_iteration):
        training_episode(tf_jungle_env, agent_list)
        a_return_list = [a.episode_return() for a in agent_list]

        games.append({
            'iteration': iteration,
            'game': game,
            'final_step': tf_jungle_env.current_time_step(),
            'a2_return': a_return_list[0],
            'a3_return': a_return_list[1],
            'a4_return': a_return_list[2],
            'a5_return': a_return_list[3]

        })


def train():
    losses = []
    for _ in range(train_steps_per_iteration):
        for a in agent_list:
            losses.append(a.train_iteration().loss.numpy())

        loss_infos.append({
            'iteration': iteration,
            'agent2': losses[0],
            'agent3': losses[1],
            'agent4': losses[2],
            'agent5': losses[3]
        })


def plot_history():

    games_data = pd.DataFrame.from_records(games)
    loss_data = pd.DataFrame.from_records(loss_infos)
    for i in range(len(agent_list)):
        loss_data[f'agent{i+2}_log'] = np.log(loss_data[f"agent{i+2}"])

    fig, axs = plt.subplots(2, 2, figsize=(15, 12))

    loss_melted = pd.melt(loss_data,
                          id_vars=['iteration'],
                          value_vars=[f'agent{i + 2}_log' for i in range(len(agent_list))])
    smoothing = iteration // 50
    loss_melted.iteration = smoothing * (loss_melted.iteration // smoothing)

    sns.lineplot(ax=axs[0][0],
                 x='iteration', hue='variable',
                 y='value', data=loss_melted)
    axs[0][0].set_title('Loss History')
    axs[0][0].set_ylabel('log-loss')

    returns_melted = pd.melt(games_data,
                             id_vars=['iteration'],
                             value_vars=['p1_return', 'p2_return'])
    returns_melted.iteration = smoothing * \
        (returns_melted.iteration // smoothing)
    sns.lineplot(ax=axs[0][1],
                 x='iteration', hue='variable',
                 y='value', data=returns_melted)
    axs[0][1].set_title('Return History')
    axs[0][1].set_ylabel('return')

    games_data['p1_win'] = games_data.outcome == 'p1_win'
    games_data['p2_win'] = games_data.outcome == 'p2_win'
    games_data['illegal'] = games_data.outcome == 'illegal'
    grouped_games_data = games_data.groupby('iteration')
    cols = ['game', 'p1_win', 'p2_win', 'illegal']
    grouped_games_data = grouped_games_data[cols]
    game_totals = grouped_games_data.max()['game'] + 1
    summed_games_data = grouped_games_data.sum()
    summed_games_data['p1_win_rate'] = summed_games_data.p1_win / game_totals
    summed_games_data['p2_win_rate'] = summed_games_data.p2_win / game_totals
    summed_games_data['illegal_rate'] = summed_games_data.illegal / game_totals
    summed_games_data['iteration'] = smoothing * \
        (summed_games_data.index // smoothing)

    sns.lineplot(ax=axs[1][0],
                 x='iteration',
                 y='p1_win_rate',
                 data=summed_games_data,
                 label='Player 1 Win Rate')
    sns.lineplot(ax=axs[1][0],
                 x='iteration',
                 y='p2_win_rate',
                 data=summed_games_data,
                 label='Player 2 Win Rate')
    sns.lineplot(ax=axs[1][0],
                 x='iteration',
                 y='illegal_rate',
                 data=summed_games_data,
                 label='Illegal Ending Rate')
    axs[1][0].set_title('Outcomes History')
    axs[1][0].set_ylabel('ratio')

    plt.show()
#


try:
    while iteration < num_iterations:
        collect_training_data()
        train()
        iteration += 1
        # if iteration % plot_interval == 0:
        #     plot_history()
        #     clear_output(wait=True)

except KeyboardInterrupt:
    # clear_output(wait=True)
    print('Interrupting training, plotting history...')
    # plot_history()
