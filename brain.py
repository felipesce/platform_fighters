from random import random
import math

from vpython import vec

def relu(x):
    return x * (0.1 if x<0 else 1.0)

def logistic(x):
    if x > 500:
        return 1
    elif x < -500:
        return 0
    else:
        return 1 / (1 + math.exp(-x))

def get_color(value, min_th=-2.0, max_th=2.0):
    # Normalize value between -1 and 1
    normalized = max(min_th, min(max_th, value)) / max_th
    if value >= 0:
        intensity = normalized  # Scale from 0 to 1
        return vec(1-intensity, 1, 1-intensity)  # White to Green
    else:
        intensity = abs(normalized)  # Scale from 0 to 1
        return vec(1, 1-intensity, 1-intensity)  # White to Red

class Neuron:

    def __init__(self, 
        activation_function = relu,
        y = 0
        
        ):

        self.F = activation_function
        self.preval = 0.0
        self.val = 0.0
        self.kind = "INPUT"

        self.y = y

        self.connections = []

    def connect(self, other, w=1.0):
        self.connections.append([other, w])

    def propagate(self):
        for c in self.connections:
            c[0].preval += c[1] * self.val

    def activate(self):
        self.val = self.F(self.preval)

class Brain:

    def __init__(self, N_in, N_out):
        
        self.N_in = N_in
        self.N_out = N_out

        self.N_tot = N_in + N_out

        # Input neurons use relu
        self.input_neurons = [Neuron(relu, n) for n in range(N_in)]
        # Output neurons use logistic function to get 0-1 range
        self.output_neurons = [Neuron(logistic, n) for n in range(N_out)]

        # Initialize weights between -1 and 1
        for i in self.input_neurons:
            for j in self.output_neurons:
                i.connect(j, random() * 2 - 1)

    def forward(self, inputs):
        # Reset all neurons
        for n in self.input_neurons + self.output_neurons:
            n.preval = 0.0
            n.val = 0.0
        
        # Set input values
        for n, inp in zip(self.input_neurons, inputs):
            n.val = inp
            n.propagate()
        
        # Activate output neurons
        for n in self.output_neurons:
            n.activate()
        
        # Return output values
        return [n.val for n in self.output_neurons]

if __name__ == "__main__":
    from vpython import *

    brain = Brain(3,4)

    # Create visuals and add labels
    input_visuals = [sphere(pos = vec(-5,2*n.y,0)) for n in brain.input_neurons]
    output_visuals = [sphere(pos = vec(5,2*n.y,0)) for n in brain.output_neurons]
    iteration_label = label(pos=vec(0,3,0), text='Iteration: 0')
    
    # Add value labels for neurons
    input_labels = [label(pos=vec(-5,2*n.y+1,0), text='0.00') for n in brain.input_neurons]
    output_labels = [label(pos=vec(5,2*n.y+1,0), text='0.00') for n in brain.output_neurons]
    
    # Create connection lines
    connection_lines = []
    for i, input_n in enumerate(brain.input_neurons):
        for j, output_n in enumerate(brain.output_neurons):
            line = cylinder(pos=input_visuals[i].pos, axis=output_visuals[j].pos-input_visuals[i].pos,
                          radius=0.05)
            connection_lines.append(line)
            # Color based on weight using the helper function
            line.color = get_color(input_n.connections[j][1])

    for i in range(100):
        rate(1)  # Slow down animation
        iteration_label.text = f'Iteration: {i}'  # Update iteration counter

        for n, visual, label in zip(brain.input_neurons, input_visuals, input_labels):
            n.val = random() * 2 - 1
            n.propagate()
            visual.color = get_color(n.val)
            label.text = f'{n.val:.2f}'

        for n, visual, label in zip(brain.output_neurons, output_visuals, output_labels):
            n.activate()
            visual.color = get_color(n.val)
            label.text = f'{n.val:.2f}'

