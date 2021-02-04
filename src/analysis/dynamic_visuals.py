from common.rotations import ned_to_body
from common.tools import angle_of_attack, angle_of_sideslip, speed
from stl import mesh
from matplotlib import animation, pyplot as plt
from mpl_toolkits import mplot3d
from numpy import array, append, cos, deg2rad, interp, linalg, linspace, pi, sin


def aircraft(t_in, cg, b, ub, vb, wb, phi_deg, theta_deg, psi_deg, nz, roll_cmd, pitch_cmd, yaw_cmd,
             roll_out, pitch_out, yaw_out, dt=0.1, alpha=0.2, kv=2, r=25):
    t = linspace(0, t_in[-1], int(t_in[-1]/dt+1))
    phi = interp(t, t_in, deg2rad(phi_deg))
    theta = interp(t, t_in, deg2rad(theta_deg))
    psi = interp(t, t_in, deg2rad(psi_deg))

    ub = interp(t, t_in, ub)
    vb = interp(t, t_in, vb)
    wb = interp(t, t_in, wb)

    ub = append([0], ub)
    vb = append([0], vb)
    wb = append([0], wb)

    roll_out = interp(t, t_in, roll_out)
    pitch_out = interp(t, t_in, pitch_out)
    yaw_out = interp(t, t_in, yaw_out)

    roll_out = append([0], roll_out)
    pitch_out = append([0], pitch_out)
    yaw_out = append([0], yaw_out)

    roll_cmd = interp(t, t_in, roll_cmd)
    pitch_cmd = interp(t, t_in, pitch_cmd)
    yaw_cmd = interp(t, t_in, yaw_cmd)

    roll_cmd = append([0], roll_cmd)
    pitch_cmd = append([0], pitch_cmd)
    yaw_cmd = append([0], yaw_cmd)

    nz = interp(t, t_in, nz)
    nz = append([0], nz)

    t = append([0], t)
    phi = append([0], phi)
    theta = append([0], theta)
    psi = append([0], psi)

    frames = len(t)-1
    fps = 1 / dt

    # First set up the figure, the axis, and the plot element we want to animate
    figure = plt.figure(figsize=(8, 7))
    axes = mplot3d.Axes3D(figure)
    axes2 = figure.add_axes([0.75, 0.75, 0.2, 0.2])
    axes.set_xlim3d(-b, b)
    axes.set_ylim3d(-b, b)
    axes.set_zlim3d(-b, b)

    axes2.set_xlim(-30, 30)
    axes2.set_ylim(-30, 30)
    axes2.grid()

    # Load the STL files and add the vectors to the plot
    your_mesh = mesh.Mesh.from_file(r"C:\Users\Andres\Documents\Airplane.stl")

    angle = linspace(0, 2 * pi, 201)
    x_horizon = r * cos(angle)
    y_horizon = r * sin(angle)
    z_horizon = angle * 0

    x_phi = 0 * angle
    y_phi = r * cos(angle)
    z_phi = r * sin(angle)

    x_theta = r * cos(angle)
    y_theta = 0 * angle
    z_theta = r * sin(angle)

    # animation function.  This is called sequentially

    def animate(index):
        j = index
        i = index + 1
        axes.clear()
        axes2.clear()

        axes.set_xlim3d(-b/2, b/2)
        axes.set_ylim3d(-b/2, b/2)
        axes.set_zlim3d(-b/2, b/2)

        axes2.set_xlim(-30, 30)
        axes2.set_ylim(-30, 30)

        phi_i = phi[i]-phi[j]
        theta_i = (theta[i] - theta[j])
        psi_i = psi[i]-psi[j]
        r_n2b = ned_to_body(phi[i], theta[i], psi[i])
        v_n = linalg.inv(r_n2b) @ array([[ub[i]], [vb[i]], [wb[i]]])
        xv = -v_n[0] * kv / b
        yv = v_n[1] * kv / b
        zv = -v_n[2] * kv / b
        axes.quiver(cg[0], cg[1], cg[2], xv, yv, zv, color='g')
        your_mesh.rotate([1, 0, 0], phi_i, point=cg)
        your_mesh.rotate([0, 1, 0], -theta_i, point=cg)
        your_mesh.rotate([0, 0, 1], psi_i, point=cg)

        obj = mplot3d.art3d.Poly3DCollection(your_mesh.vectors)
        obj.set_alpha(alpha)
        obj.set_facecolor('b')
        axes.add_collection3d(obj)

        # Auto scale to the mesh size
        scale = your_mesh.points.flatten()
        axes.auto_scale_xyz(scale, scale, scale)
        axes.scatter3D(cg[0], cg[1], cg[2], marker='o', edgecolors='r')

        axes.scatter3D(-r+cg[0], cg[1], cg[2], marker='o', edgecolors='brown')
        axes.scatter3D(r + cg[0], cg[1], 0 + cg[2], marker='o', edgecolors='brown')
        axes.scatter3D(cg[0], r + cg[1], 0 + cg[2], marker='o', edgecolors='brown')
        axes.scatter3D(cg[0], -r + cg[1], 0 + cg[2], marker='o', edgecolors='brown')
        axes.text3D(-r+cg[0]-5, cg[1], 0+cg[2], 'N')
        axes.text3D(cg[0], r + cg[1] + 5, 0 + cg[2], 'E')
        axes.text3D(cg[0], -r + cg[1] - 5, 0 + cg[2], 'W')
        axes.text3D(r + cg[0] + 5, cg[1], 0 + cg[2], 'S')

        v_fps = speed(array([ub[i], vb[i], wb[i]]))
        aoa = angle_of_attack(array([ub[i], vb[i], wb[i]]))
        aos = angle_of_sideslip(array([ub[i], vb[i], wb[i]]))

        r_phi = r_n2b @ array([x_phi, y_phi, z_phi])
        r_theta = r_n2b @ array([x_theta, y_theta, z_theta])

        axes.plot(x_horizon+cg[0], y_horizon+cg[1], z_horizon+cg[2], color='brown')
        axes.plot(r_phi[0, :] + cg[0], r_phi[1, :] + cg[1], -r_phi[2, :] + cg[2], color='red')
        axes.plot(r_theta[0, :] + cg[0], r_theta[1, :] + cg[1], r_theta[2, :] + cg[2], color='blue')

        axes.text2D(0.05, 0.85, "V %4.1f fps \n\u03B1 %2.2f\u00B0\n \u03B2 %2.2f\u00B0\n Nz %2.2f g \n Vy %3.2f fpm"
                    % (v_fps, aoa, aos, nz[i], -v_n[2]),
                    transform=axes.transAxes)

        axes2.scatter(roll_cmd[i], pitch_cmd[i], marker='o', edgecolors='r', facecolor=None)
        axes2.scatter(yaw_cmd[i], 0, marker='^', edgecolors='b', facecolor=None)
        axes2.scatter(roll_out[i], pitch_out[i], marker='o', edgecolors='r', facecolor='k')
        axes2.scatter(yaw_out[i], 0, marker='^', edgecolors='b', facecolor='k')

        return axes, axes2

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(figure, animate, frames=frames, interval=dt/fps*1000, repeat=False,
                                   save_count=frames)
    anim.save('basic_animation.mp4', fps=int(fps), extra_args=['-vcodec', 'libx264'])
