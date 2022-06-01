
import tensorflow as tf
import numpy as np

from tf_agents.environments import py_environment
# from tf_agents.environments import tf_environment
# from tf_agents.environments import tf_py_environment
# from tf_agents.environments import utils
from tf_agents.specs import array_spec
from tf_agents.trajectories import time_step as ts
from copy import copy


class JungleEnv(py_environment.PyEnvironment):

    def __init__(self, grid=(20, 20), L=8, N=4, MAX_STEP=200):
        super(JungleEnv, self).__init__(handle_auto_reset=True)
        self.ILLEGAL_MOVE_REWARD = np.asarray(-10.0, dtype=np.float32)
        self.GETTING_FOOD_REWARD = np.asarray(1, dtype=np.float32)
        self.ATTACK_ENEMY_REWARD = np.asarray(3, dtype=np.float32)
        self.GETTING_ATTACKED_REWARD = np.asarray(-4, dtype=np.float32)
        self.FAILED_ATTACK_REWARD = np.asarray(-1.0, dtype=np.float32)
        self.NULL_REWARD = np.asarray(0.0, dtype=np.float32)

        super(JungleEnv, self).__init__(handle_auto_reset=False)
        self.grid = grid
        self.L = L
        self.N = N
        self.MAX_STEP = MAX_STEP
        self.step_counter = 0
        # actions are either [go, attack] --> [west, east, south, north]
        # i.e. 0 -> go west, 1 -> go east, ..., 4 -> attack west, 5 -> attack east, ...
        self._action_spec = array_spec.BoundedArraySpec(
            shape=(), dtype=np.int32, minimum=0, maximum=7, name='action')

        # observation is a 3x3 grid whith 0 -> nothing, 1 -> food, 2+ -> another agent
        # self._observation_spec = array_spec.BoundedArraySpec(
        #     shape=grid, dtype=np.int32, minimum=0, maximum=self.N+2, name='observation')
        self._observation_spec = array_spec.BoundedArraySpec(
            shape=(3, 3), dtype=np.int32, minimum=0, maximum=self.N+2, name='observation')

        # agent spec for separating the agents
        self._agent_spec = array_spec.BoundedArraySpec(
            (), np.int32, minimum=2, maximum=self.N + 2, name='agent')

        # we initialize the state with L foods and N agents randomly in thw grid
        self._episode_ended = False
        self._reset()

    def action_spec(self):
        return {
            'action': self._action_spec,
            'agent': self._agent_spec
        }

    def observation_spec(self):
        return self._observation_spec

    def _reset(self):
        self.step_counter = 0

        # we initialize the state with L foods and N agents randomly in the grid
        self._state = np.zeros(shape=self.grid, dtype=np.int32)
        reward_record = np.zeros((self.grid[1], ), np.int32)
        r1 = np.random.permutation(range(self.grid[0]))[:self.L]
        r2 = np.random.permutation(range(self.grid[1]))[:self.L]
        self.food_locs = np.c_[r1, r2]
        self._state[self.food_locs[:, 0], self.food_locs[:, 1]] = 1
        # place agents where there is no food
        free_cells = np.where(self._state == 0)
        r = np.random.permutation(range(free_cells[0].shape[0]))[:self.N]
        self.agent_locs = np.c_[free_cells[0][r], free_cells[1][r]]
        # each agent has a unique id between 2 and N + 2
        for i in range(self.N):
            self._state[self.agent_locs[i, 0], self.agent_locs[i, 1]] = i + 2

        self._state = np.c_[self._state, reward_record]
        self.update_periodic_state()
        self.reset_attack_records()
        #
        self._episode_ended = False
        return ts.restart(self._state)

    def reset_attack_records(self):
        self.attack_reports = [False] * self.N

    def _step(self, action):
        self.step_counter += 1
        if self.step_counter >= self.MAX_STEP:
            self._episode_ended = True
            self._reset()

            return ts.termination(self._state, reward=0.0)

        possible_actions = [self.move_west, self.move_east, self.move_south, self.move_north,
                            self.attack_west, self.attack_east, self.attack_south, self.attack_north]
        if self._episode_ended:
            # The last action ended the episode. Ignore the current action and start
            # a new episode.
            return self.reset()

        agent_id = action["agent"]
        reward = 0.0
        if self.attack_reports[agent_id - 2]:
            reward += self.GETTING_ATTACKED_REWARD

        if (self.step_counter % self.N) == 0:
            self.reset_attack_records()
        agent_loc = self.agent_locs[agent_id - 2, :]
        r, is_illegal = possible_actions[action["action"]](
            agent_loc, agent_id)
        reward += r
        if is_illegal:
            self._reset()
            self._episode_ended = True
            self._state[agent_id - 2, -1] += int(self.ILLEGAL_MOVE_REWARD)
            return ts.termination(self._state, reward=copy(reward))
        else:
            self._state[agent_id - 2, -1] += int(reward)

            return ts.transition(self._state, reward=copy(reward), discount=1.0)

    def check_food(self, loc):
        l = [
            self._periodic_state[loc[0] + 1, loc[1]],
            self._periodic_state[loc[0] - 1, loc[1]],
            self._periodic_state[loc[0], loc[1] + 1],
            self._periodic_state[loc[0], loc[1] - 1],
        ]
        if 1 in l:
            return True
        else:
            return False

    def update_periodic_state(self):
        _ = np.r_[self._state[:, :-1], self._state[:, :-1]]
        self._periodic_state = np.c_[_, _]

    def move_west(self, pos, id):
        x = pos[0]
        y = pos[1]
        if self._periodic_state[x - 1, y] == 0:
            self._state[:, :-1][pos[0], pos[1]] = 0
            new_x = pos[0] - 1 if pos[0] != 0 else self.grid[0] - 1
            self._state[:, :-1][new_x, pos[1]] = id
            self.update_periodic_state()
            self.agent_locs[id - 2, 0] = new_x
            # check for food near by
            x = new_x
            food = self.check_food([x, y])
            if food:
                return self.GETTING_FOOD_REWARD, False
            else:
                return self.NULL_REWARD, False

        else:
            return self.ILLEGAL_MOVE_REWARD, True

    def move_east(self, pos, id):
        x = pos[0]
        y = pos[1]
        if self._periodic_state[x + 1, y] == 0:
            self._state[:, :-1][pos[0], pos[1]] = 0
            new_x = pos[0] + 1 if pos[0] != self.grid[0] - 1 else 0
            self._state[:, :-1][new_x, pos[1]] = id
            self.update_periodic_state()
            self.agent_locs[id - 2, 0] = new_x
            # check for food near by
            x = new_x
            food = self.check_food([x, y])
            if food:
                return self.GETTING_FOOD_REWARD, False
            else:
                return self.NULL_REWARD, False

        else:
            return self.ILLEGAL_MOVE_REWARD, True

    def move_south(self, pos, id):
        x = pos[0]
        y = pos[1]
        if self._periodic_state[x, y - 1] == 0:
            self._state[:, :-1][pos[0], pos[1]] = 0
            new_y = pos[1] - 1 if pos[1] != 0 else self.grid[1] - 1
            self._state[:, :-1][pos[0], new_y] = id
            self.update_periodic_state()
            self.agent_locs[id - 2, 1] = new_y
            # check for food near by
            y = new_y
            food = self.check_food([x, y])
            if food:
                return self.GETTING_FOOD_REWARD, False
            else:
                return self.NULL_REWARD, False
        else:
            return self.ILLEGAL_MOVE_REWARD, True

    def move_north(self, pos, id):
        x = pos[0]
        y = pos[1]
        if self._periodic_state[x, y + 1] == 0:
            self._state[:, :-1][pos[0], pos[1]] = 0
            new_y = pos[1] + 1 if pos[1] != self.grid[1] - 1 else 0
            self._state[:, :-1][pos[0], new_y] = id
            self.update_periodic_state()
            self.agent_locs[id - 2, 1] = new_y
            # check for food near by
            y = new_y
            food = self.check_food([x, y])
            if food:
                return self.GETTING_FOOD_REWARD, False
            else:
                return self.NULL_REWARD, False
        else:
            return self.ILLEGAL_MOVE_REWARD, True

    def attack_west(self, pos, id):
        x = pos[0]
        y = pos[1]
        if (self._periodic_state[x - 1, y] != 0) and (self._periodic_state[x - 1, y] != 1):
            attacked_enemy_id = self._periodic_state[x - 1, y]
            self.attack_reports[attacked_enemy_id - 2] = True
            return self.ATTACK_ENEMY_REWARD, False
        else:
            return self.FAILED_ATTACK_REWARD, False

    def attack_east(self, pos, id):
        x = pos[0]
        y = pos[1]
        if (self._periodic_state[x + 1, y] != 0) and (self._periodic_state[x + 1, y] != 1):
            attacked_enemy_id = self._periodic_state[x + 1, y]
            self.attack_reports[attacked_enemy_id - 2] = True
            return self.ATTACK_ENEMY_REWARD, False
        else:
            return self.FAILED_ATTACK_REWARD, False

    def attack_south(self, pos, id):
        x = pos[0]
        y = pos[1]
        if (self._periodic_state[x, y - 1] != 0) and (self._periodic_state[x, y - 1] != 1):
            attacked_enemy_id = self._periodic_state[x, y - 1]
            self.attack_reports[attacked_enemy_id - 2] = True
            return self.ATTACK_ENEMY_REWARD, False
        else:
            return self.FAILED_ATTACK_REWARD, False

    def attack_north(self, pos, id):
        x = pos[0]
        y = pos[1]
        if (self._periodic_state[x, y + 1] != 0) and (self._periodic_state[x, y + 1] != 1):
            attacked_enemy_id = self._periodic_state[x, y + 1]
            self.attack_reports[attacked_enemy_id - 2] = True
            return self.ATTACK_ENEMY_REWARD, False
        else:
            return self.FAILED_ATTACK_REWARD, False
