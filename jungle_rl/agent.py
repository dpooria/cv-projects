
import numpy as np
import tensorflow as tf

from tf_agents.agents import DqnAgent
from tf_agents.agents.tf_agent import LossInfo
from tf_agents.environments import TFPyEnvironment
from tf_agents.specs import TensorSpec
from typing import Callable
from tf_agents.trajectories.time_step import TimeStep
from tf_agents.replay_buffers.tf_uniform_replay_buffer import TFUniformReplayBuffer
from tf_agents.utils import common
from tf_agents.networks.q_network import QNetwork
from tf_agents.trajectories import trajectory
from tf_agents.trajectories.trajectory import Trajectory

from jenv import JungleEnv


class IMAgent(DqnAgent):

    def __init__(self,
                 env: JungleEnv,
                 observation_spec: TensorSpec = None,
                 action_spec: TensorSpec = None,
                 action_fn: Callable = lambda action: action,
                 name: str = 'IMAgent',
                 id: int = 2,
                 q_network=None,
                 # training params
                 replay_buffer_max_length: int = 1000,
                 learning_rate: float = 1e-5,
                 training_batch_size: int = 8,
                 training_parallel_calls: int = 3,
                 training_prefetch_buffer_size: int = 3,
                 training_num_steps: int = 2,
                 **dqn_kwargs):

        self._env = env
        # self._reward_fn = reward_fn
        self._name = name
        self._observation_spec = observation_spec or self._env.observation_spec()
        self._action_spec = action_spec or self._env.action_spec()
        self._action_fn = action_fn
        self.id = id
        q_network = q_network or self._build_q_net()
        env_ts_spec = self._env.time_step_spec()
        time_step_spec = TimeStep(
            step_type=env_ts_spec.step_type,
            reward=env_ts_spec.reward,
            discount=env_ts_spec.discount,
            observation=q_network.input_tensor_spec
        )

        optimizer = tf.keras.optimizers.Adam(learning_rate=learning_rate)
        super().__init__(time_step_spec,
                         self._action_spec,
                         q_network,
                         optimizer,
                         name=name,
                         **dqn_kwargs)

        self._policy_state = self.policy.get_initial_state(
            batch_size=self._env.batch_size)
        self._rewards = []

        self._replay_buffer = TFUniformReplayBuffer(
            data_spec=self.collect_data_spec,
            batch_size=self._env.batch_size,
            max_length=replay_buffer_max_length)

        self._training_batch_size = training_batch_size
        self._training_parallel_calls = training_parallel_calls
        self._training_prefetch_buffer_size = training_prefetch_buffer_size
        self._training_num_steps = training_num_steps
        self.train = common.function(self.train)

    def _build_q_net(self):
        fc_layer_params = (50,)
        q_net = QNetwork(
            self._observation_spec,
            self._action_spec,
            fc_layer_params=fc_layer_params)

        q_net.create_variables()
        q_net.summary()
        return q_net

    def reset(self):
        self._policy_state = self.policy.get_initial_state(
            batch_size=self._env.batch_size
        )
        self._rewards = []

    def episode_return(self) -> float:
        return np.sum(self._rewards)

    def _reward_fn(self, time_step):
        observation = time_step.observation
        self.reward = observation[:, self.id - 2, -1]
        return tf.cast(self.reward, tf.float32)

    def _observation_fn(self, observation: tf.Tensor) -> tf.Tensor:
        # each agent only sees a 3x3 grid according to it's current position
        agent_loc = list(np.where(observation[:, :, :-1] == self.id))
        agent_loc = [agent_loc[0][0], agent_loc[1][0], agent_loc[2][0]]
        self.grid = tf.cast(tf.shape(observation[0, :, :-1]), dtype=tf.int64)

        state_periodic = observation[:, :, :-1]

        state_periodic = tf.concat((state_periodic, state_periodic), axis=1)
        state_periodic = tf.concat([state_periodic, state_periodic], axis=2)
        state_periodic = tf.reshape(
            state_periodic, [-1, state_periodic.shape[1], state_periodic.shape[2]])
        # state_periodic = self._env._periodic_step
        loc_periodic = agent_loc
        a = 0
        b = 0
        if agent_loc[1] == 0:
            a = agent_loc[1] + self.grid[0]
        if agent_loc[2] == 0:
            b = agent_loc[2] + self.grid[1]

        loc_periodic[1] += a
        loc_periodic[2] += b
        observation = state_periodic[:, loc_periodic[1] -
                                     1:loc_periodic[1] + 2, loc_periodic[2] - 1:loc_periodic[2] + 2]

        return observation

    def _augment_time_step(self, time_step: TimeStep) -> TimeStep:
        reward = self._reward_fn(time_step)
        reward = tf.convert_to_tensor(reward, dtype=np.float32)
        if reward.shape != time_step.reward.shape:
            reward = tf.reshape(reward, time_step.reward.shape)

        observation = self._observation_fn(time_step.observation)

        return TimeStep(
            step_type=time_step.step_type,
            reward=reward,
            discount=time_step.discount,
            observation=observation
        )

    def _current_time_step(self) -> TimeStep:
        time_step = self._env.current_time_step()
        time_step = self._augment_time_step(time_step)
        return time_step

    def _step_environment(self, action) -> TimeStep:
        action = self._action_fn(action)
        time_step = self._env.step(action)
        time_step = self._augment_time_step(time_step)
        return time_step

    def act(self, collect=False) -> Trajectory:
        time_step = self._current_time_step()
        if collect:
            policy_step = self.collect_policy.action(
                time_step, policy_state=self._policy_state)
        else:
            policy_step = self.policy.action(
                time_step, policy_state=self._policy_state)

        self._policy_state = policy_step.state
        next_time_step = self._step_environment(policy_step.action)
        traj = trajectory.from_transition(
            time_step, policy_step, next_time_step)

        self._rewards.append(next_time_step.reward)

        if collect:
            self._replay_buffer.add_batch(traj)

        return traj

    def train_iteration(self) -> LossInfo:
        experience, info = self._replay_buffer.get_next(
            sample_batch_size=self._training_batch_size,
            num_steps=self._training_num_steps
        )
        return self.train(experience)
