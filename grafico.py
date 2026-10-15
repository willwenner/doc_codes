import numpy as np
import plotly.graph_objects as go

R, r = 3.0, 1.0

theta = np.linspace(0, 2 * np.pi, 60)
phi = np.linspace(0, 2 * np.pi, 60)
THETA, PHI = np.meshgrid(theta, phi)

X = (R + r * np.cos(THETA)) * np.cos(PHI)
Y = (R + r * np.cos(THETA)) * np.sin(PHI)
Z = r * np.sin(THETA)

theta_pol = np.linspace(0, 2 * np.pi, 100)
x_pol = R + r * np.cos(theta_pol)
y_pol = np.zeros_like(theta_pol)
z_pol = r * np.sin(theta_pol)

fig = go.Figure()

fig.add_trace(go.Surface(x=X, y=Y, z=Z, opacity=0.5, colorscale='Blues', showscale=False))
fig.add_trace(go.Scatter3d(x=x_pol, y=y_pol, z=z_pol, mode='lines', line=dict(color='red', width=8)))

axis_config = dict(
    showgrid=False, zeroline=True, showline=True, ticks='', showticklabels=False, showbackground=False
)

# --- DEFININDO LIMITES DE ZOOM E POSIÇÃO DA CÂMERA ---
fig.update_layout(
    scene=dict(
        xaxis=dict(axis_config, range=[1.8, 4.2]),
        yaxis=dict(axis_config, range=[-1.5, 1.5]),
        zaxis=dict(axis_config, range=[-1.5, 1.5]),
        aspectmode='cube',
        camera=dict(
            eye=dict(x=1.2, y=1.2, z=0.8),  # Posição de aproximação da câmera
            center=dict(x=0.3, y=0, z=0)     # Aponta a câmera para a seção
        )
    )
)

fig.show()