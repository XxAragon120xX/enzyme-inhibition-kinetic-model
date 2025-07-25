import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import tkinter as tk
from tkinter import ttk, messagebox
import seaborn as sns
from matplotlib.gridspec import GridSpec


plt.style.use('default')  
sns.set_palette("husl")

class InteractiveKineticModel:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("Modelo Cinético de Inhibición Enzimática - Interactivo")
        self.root.geometry("1400x900")
        self.root.configure(bg='#f0f0f0')
        
        # Parámetros del modelo
        self.parameters = {
            'mu_max': tk.DoubleVar(value=0.95),
            'Ki_comp1': tk.DoubleVar(value=0.12),
            'Ki_comp2': tk.DoubleVar(value=10.0),
            'K_max': tk.DoubleVar(value=1e9),
            'ka': tk.DoubleVar(value=0.1),
            'ke': tk.DoubleVar(value=0.05),
            'initial_conc': tk.DoubleVar(value=0.2),
            'time_points': tk.IntVar(value=100)
        }
        
        # Variable para seleccionar compuesto
        self.selected_compound = tk.StringVar(value="compound1")
        
        # Datos experimentales
        self.experimental_data = {
            'compound1': pd.DataFrame([
                {'conc': 0.00, 'mu': 0.95, 'g': 0.73, 'viable': 1.0},
                {'conc': 0.05, 'mu': 0.85, 'g': 0.82, 'viable': 0.89},
                {'conc': 0.10, 'mu': 0.75, 'g': 0.92, 'viable': 0.79},
                {'conc': 0.15, 'mu': 0.55, 'g': 1.26, 'viable': 0.58},
                {'conc': 0.20, 'mu': 0.40, 'g': 1.73, 'viable': 0.42}
            ]),
            'compound2': pd.DataFrame([
                {'conc': 0.00, 'mu': 0.95, 'g': 0.73, 'viable': 1.0},
                {'conc': 0.05, 'mu': 0.90, 'g': 0.77, 'viable': 0.95},
                {'conc': 0.10, 'mu': 0.91, 'g': 0.76, 'viable': 0.96},
                {'conc': 0.15, 'mu': 0.89, 'g': 0.78, 'viable': 0.94},
                {'conc': 0.20, 'mu': 0.90, 'g': 0.77, 'viable': 0.95}
            ])
        }
        
        self.setup_ui()
        self.bind_events()
        self.update_plots()
        
    def setup_ui(self):
        """Configurar la interfaz de usuario"""
        # Frame principal
        main_frame = ttk.Frame(self.root)
        main_frame.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)
        
        # Configurar grid
        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_columnconfigure(0, weight=1)
        main_frame.grid_rowconfigure(1, weight=1)
        main_frame.grid_columnconfigure(1, weight=1)
        
        # Título
        title_frame = ttk.Frame(main_frame)
        title_frame.grid(row=0, column=0, columnspan=2, sticky="ew", pady=(0, 10))
        
        title_label = ttk.Label(title_frame, 
                               text="Modelo Cinético de Inhibición Enzimática", 
                               style="Title.TLabel")
        title_label.pack()
        
        subtitle_label = ttk.Label(title_frame,
                                 text="Simulación interactiva del efecto de inhibidores de DHFR sobre el crecimiento de E. coli",
                                 style="Subtitle.TLabel")
        subtitle_label.pack()
        
        # Panel de controles
        self.setup_control_panel(main_frame)
        
        # Panel de gráficos
        self.setup_plot_panel(main_frame)
        
        # Panel de resultados
        self.setup_results_panel(main_frame)
        
        # Configurar estilos
        self.setup_styles()
    
    def setup_styles(self):
        """Configurar estilos personalizados"""
        style = ttk.Style()
        
        # Estilo para título
        style.configure("Title.TLabel", 
                       font=("Arial", 18, "bold"),
                       foreground="#2c3e50")
        
        style.configure("Subtitle.TLabel",
                       font=("Arial", 12),
                       foreground="#7f8c8d")
        
        # Estilo para secciones
        style.configure("Section.TLabel",
                       font=("Arial", 14, "bold"),
                       foreground="#34495e")
        
        # Estilo para parámetros
        style.configure("Param.TLabel",
                       font=("Arial", 10),
                       foreground="#2c3e50")
        
    def setup_control_panel(self, parent):
        """Configurar panel de controles"""
        control_frame = ttk.LabelFrame(parent, text="Parámetros del Modelo", padding=10)
        control_frame.grid(row=1, column=0, sticky="nsew", padx=(0, 5))
        
        # Selector de compuesto
        compound_frame = ttk.LabelFrame(control_frame, text="Seleccionar Compuesto", padding=10)
        compound_frame.pack(fill="x", pady=(0, 10))
        
        ttk.Radiobutton(compound_frame, text="Compuesto 1 (Alta actividad)", 
                       variable=self.selected_compound, value="compound1",
                       command=self.on_compound_change).pack(anchor="w")
        
        ttk.Radiobutton(compound_frame, text="Compuesto 2 (Baja actividad)", 
                       variable=self.selected_compound, value="compound2",
                       command=self.on_compound_change).pack(anchor="w")
        
        # Parámetros de inhibición
        inhibition_frame = ttk.LabelFrame(control_frame, text="Parámetros de Inhibición", padding=10)
        inhibition_frame.pack(fill="x", pady=(0, 10))
        
        self.create_parameter_control(inhibition_frame, "μ máx (h⁻¹):", 
                                    self.parameters['mu_max'], 0.1, 2.0, 0.01)
        
        self.create_parameter_control(inhibition_frame, "Ki Comp1 (mmol/L):", 
                                    self.parameters['Ki_comp1'], 0.01, 1.0, 0.01)
        
        self.create_parameter_control(inhibition_frame, "Ki Comp2 (mmol/L):", 
                                    self.parameters['Ki_comp2'], 1.0, 50.0, 0.5)
        
        # Parámetros farmacocinéticos
        pharmaco_frame = ttk.LabelFrame(control_frame, text="Parámetros Farmacocinéticos", padding=10)
        pharmaco_frame.pack(fill="x", pady=(0, 10))
        
        self.create_parameter_control(pharmaco_frame, "ka (h⁻¹):", 
                                    self.parameters['ka'], 0.01, 0.5, 0.01)
        
        self.create_parameter_control(pharmaco_frame, "ke (h⁻¹):", 
                                    self.parameters['ke'], 0.01, 0.2, 0.01)
        
        self.create_parameter_control(pharmaco_frame, "Conc. inicial (mmol/L):", 
                                    self.parameters['initial_conc'], 0.05, 0.5, 0.01)
        
        # Parámetros de simulación
        sim_frame = ttk.LabelFrame(control_frame, text="Parámetros de Simulación", padding=10)
        sim_frame.pack(fill="x", pady=(0, 10))
        
        self.create_parameter_control(sim_frame, "Puntos de tiempo:", 
                                    self.parameters['time_points'], 50, 200, 10, is_int=True)
        
        # Botones de acción
        button_frame = ttk.Frame(control_frame)
        button_frame.pack(fill="x", pady=10)
        
        ttk.Button(button_frame, text="Actualizar Gráficos", 
                  command=self.update_plots).pack(side="left", padx=(0, 5))
        
        ttk.Button(button_frame, text="Resetear Parámetros", 
                  command=self.reset_parameters).pack(side="left", padx=5)
        
        ttk.Button(button_frame, text="Exportar Datos", 
                  command=self.export_data).pack(side="left", padx=5)
        
    def create_parameter_control(self, parent, label, variable, min_val, max_val, step, is_int=False):
        """Crear control para parámetro"""
        frame = ttk.Frame(parent)
        frame.pack(fill="x", pady=2)
        
        ttk.Label(frame, text=label, style="Param.TLabel").pack(side="left")
        
        if is_int:
            scale = ttk.Scale(frame, from_=min_val, to=max_val, 
                            variable=variable, orient="horizontal",
                            command=lambda x: self.on_parameter_change())
        else:
            scale = ttk.Scale(frame, from_=min_val, to=max_val, 
                            variable=variable, orient="horizontal",
                            command=lambda x: self.on_parameter_change())
        
        scale.pack(side="left", fill="x", expand=True, padx=5)
        
        value_label = ttk.Label(frame, text=f"{variable.get():.3f}")
        value_label.pack(side="right")
        
        # Actualizar etiqueta cuando cambie el valor
        def update_label(*args):
            if is_int:
                value_label.config(text=f"{int(variable.get())}")
            else:
                value_label.config(text=f"{variable.get():.3f}")
        
        variable.trace("w", update_label)
        
    def setup_plot_panel(self, parent):
        """Configurar panel de gráficos"""
        plot_frame = ttk.LabelFrame(parent, text="Visualizaciones", padding=5)
        plot_frame.grid(row=1, column=1, sticky="nsew", padx=(5, 0))
        
        # Notebook para múltiples pestañas
        self.notebook = ttk.Notebook(plot_frame)
        self.notebook.pack(fill="both", expand=True)
        
        # Pestaña de inhibición
        self.inhibition_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.inhibition_frame, text="Curvas de Inhibición")
        
        # Pestaña de cinética
        self.kinetic_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.kinetic_frame, text="Cinética Integrada")
        
        # Pestaña de comparación
        self.comparison_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.comparison_frame, text="Análisis Comparativo")
        
        self.setup_plots()
        
    def setup_plots(self):
        """Configurar las figuras de matplotlib"""
        # Figura para curvas de inhibición
        self.fig_inhibition = Figure(figsize=(10, 6), dpi=100)
        self.canvas_inhibition = FigureCanvasTkAgg(self.fig_inhibition, self.inhibition_frame)
        self.canvas_inhibition.get_tk_widget().pack(fill="both", expand=True)
        
        # Figura para cinética integrada
        self.fig_kinetic = Figure(figsize=(10, 8), dpi=100)
        self.canvas_kinetic = FigureCanvasTkAgg(self.fig_kinetic, self.kinetic_frame)
        self.canvas_kinetic.get_tk_widget().pack(fill="both", expand=True)
        
        # Figura para comparación
        self.fig_comparison = Figure(figsize=(10, 6), dpi=100)
        self.canvas_comparison = FigureCanvasTkAgg(self.fig_comparison, self.comparison_frame)
        self.canvas_comparison.get_tk_widget().pack(fill="both", expand=True)
        
    def setup_results_panel(self, parent):
        """Configurar panel de resultados"""
        results_frame = ttk.LabelFrame(parent, text="Resultados y Análisis", padding=10)
        results_frame.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(10, 0))
        
        # Frame para resultados numéricos
        numeric_frame = ttk.Frame(results_frame)
        numeric_frame.pack(fill="x")
        
        # Variables para mostrar resultados
        self.results_vars = {
            'ed50_comp1': tk.StringVar(),
            'ed50_comp2': tk.StringVar(),
            'potencia_rel': tk.StringVar(),
            'ki_ratio': tk.StringVar()
        }
        
        # Crear etiquetas para resultados
        for i, (key, var) in enumerate(self.results_vars.items()):
            ttk.Label(numeric_frame, textvariable=var, 
                     style="Param.TLabel").grid(row=0, column=i, padx=10, sticky="w")
        
    def inhibition_model(self, concentration, mu_max, Ki):
        """Modelo de inhibición competitiva"""
        return mu_max / (1 + concentration / Ki)
    
    def pharmacokinetic_model(self, time, initial_conc, ka, ke):
        """Modelo farmacocinético"""
        k_total = ka + ke
        return initial_conc * np.exp(-k_total * time)
    
    def population_growth_model(self, time, mu, K_max, N0=1e6):
        """Modelo de crecimiento poblacional"""
        if mu <= 0:
            return N0 * np.exp(-0.1 * time)
        
        biomass = (K_max * N0) / (N0 + (K_max - N0) * np.exp(-mu * time))
        return biomass
    
    def calculate_kinetic_data(self):
        """Calcular datos cinéticos integrados"""
        dt = 8 / self.parameters['time_points'].get()
        time_array = np.arange(0, 8 + dt, dt)
        
        results = []
        for time in time_array:
            current_conc = self.pharmacokinetic_model(
                time, 
                self.parameters['initial_conc'].get(),
                self.parameters['ka'].get(),
                self.parameters['ke'].get()
            )
            
            # Usar Ki según el compuesto seleccionado
            ki_value = (self.parameters['Ki_comp1'].get() if 
                       self.selected_compound.get() == "compound1" else 
                       self.parameters['Ki_comp2'].get())
            
            current_mu = self.inhibition_model(
                current_conc,
                self.parameters['mu_max'].get(),
                ki_value
            )
            
            biomass = self.population_growth_model(
                time, current_mu, self.parameters['K_max'].get()
            )
            log_biomass = np.log10(biomass)
            
            results.append({
                'time': time,
                'concentration': current_conc,
                'mu': current_mu,
                'biomass': log_biomass,
                'raw_biomass': biomass
            })
        
        return pd.DataFrame(results)
    
    def update_inhibition_plot(self):
        """Actualizar gráfico de inhibición"""
        self.fig_inhibition.clear()
        # CORRECCIÓN: Remover el parámetro figsize
        ax1, ax2 = self.fig_inhibition.subplots(1, 2)
        
        compounds = ['compound1', 'compound2']
        titles = ['Compuesto 1 - Alta Actividad', 'Compuesto 2 - Baja Actividad']
        colors = ['red', 'green']
        ki_values = [self.parameters['Ki_comp1'].get(), self.parameters['Ki_comp2'].get()]
        
        for i, (compound, title, color, ki) in enumerate(zip(compounds, titles, colors, ki_values)):
            ax = ax1 if i == 0 else ax2
            
            # Datos experimentales
            exp_data = self.experimental_data[compound]
            ax.scatter(exp_data['conc'], exp_data['mu'], 
                      color=color, alpha=0.7, s=80, label='Datos Experimentales', zorder=5)
            
            # Curva teórica
            conc_range = np.linspace(0, 0.25, 100)
            predicted_mu = [self.inhibition_model(conc, self.parameters['mu_max'].get(), ki) 
                           for conc in conc_range]
            
            ax.plot(conc_range, predicted_mu, color='blue', linewidth=2, label='Modelo Ajustado')
            
            ax.set_xlabel('Concentración (mmol/L)')
            ax.set_ylabel('μ (h⁻¹)')
            ax.set_title(title)
            ax.legend()
            ax.grid(True, alpha=0.3)
            
            # Información de parámetros
            ed50 = ki  # Para inhibición competitiva
            ax.text(0.05, 0.95, f'Ki: {ki:.3f} mmol/L\nED50: {ed50:.3f} mmol/L',
                   transform=ax.transAxes, fontsize=9, verticalalignment='top',
                   bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.8))
        
        self.fig_inhibition.tight_layout()
        self.canvas_inhibition.draw()
    
    def update_kinetic_plot(self):
        """Actualizar gráfico de cinética integrada"""
        kinetic_data = self.calculate_kinetic_data()
        
        self.fig_kinetic.clear()
        gs = GridSpec(2, 2, figure=self.fig_kinetic, hspace=0.3, wspace=0.3)
        
        # Farmacocinética
        ax1 = self.fig_kinetic.add_subplot(gs[0, 0])
        ax1.plot(kinetic_data['time'], kinetic_data['concentration'], 
                'r-', linewidth=2, label='[Inhibidor]')
        ax1.set_xlabel('Tiempo (h)')
        ax1.set_ylabel('Concentración (mmol/L)')
        ax1.set_title('Farmacocinética del Inhibidor')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # Tasa de crecimiento
        ax2 = self.fig_kinetic.add_subplot(gs[0, 1])
        ax2.plot(kinetic_data['time'], kinetic_data['mu'], 
                'b-', linewidth=2, label='μ')
        ax2.set_xlabel('Tiempo (h)')
        ax2.set_ylabel('μ (h⁻¹)')
        ax2.set_title('Tasa de Crecimiento Específico')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        # Crecimiento de biomasa - CORRECCIÓN AQUÍ
        ax3 = self.fig_kinetic.add_subplot(gs[1, :])
        ax3.plot(kinetic_data['time'], kinetic_data['biomass'], 
                'g-', linewidth=2, label='Biomasa')
        ax3.set_xlabel('Tiempo (h)')
        ax3.set_ylabel('log₁₀(Células/mL)')
        ax3.set_title('Crecimiento de Biomasa')
        ax3.grid(True, alpha=0.3)
        ax3.legend()
        
        self.canvas_kinetic.draw()
    
    def update_comparison_plot(self):
        """Actualizar gráfico de comparación"""
        self.fig_comparison.clear()
        ax = self.fig_comparison.add_subplot(111)
        
        # Comparar ambos compuestos
        conc_range = np.linspace(0, 0.25, 100)
        
        mu_comp1 = [self.inhibition_model(conc, self.parameters['mu_max'].get(), 
                                         self.parameters['Ki_comp1'].get()) 
                   for conc in conc_range]
        
        mu_comp2 = [self.inhibition_model(conc, self.parameters['mu_max'].get(), 
                                         self.parameters['Ki_comp2'].get()) 
                   for conc in conc_range]
        
        ax.plot(conc_range, mu_comp1, 'r-', linewidth=3, label='Compuesto 1', alpha=0.8)
        ax.plot(conc_range, mu_comp2, 'g-', linewidth=3, label='Compuesto 2', alpha=0.8)
        
        # Datos experimentales
        for compound, color in [('compound1', 'red'), ('compound2', 'green')]:
            exp_data = self.experimental_data[compound]
            ax.scatter(exp_data['conc'], exp_data['mu'], 
                      color=color, alpha=0.6, s=60, zorder=5)
        
        ax.set_xlabel('Concentración (mmol/L)')
        ax.set_ylabel('μ (h⁻¹)')
        ax.set_title('Comparación de Actividad de Compuestos')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Línea horizontal para 50% de inhibición
        ax.axhline(y=self.parameters['mu_max'].get()*0.5, color='black', 
                  linestyle='--', alpha=0.5, label='50% inhibición')
        
        self.canvas_comparison.draw()
    
    def update_results(self):
        """Actualizar panel de resultados"""
        ki1 = self.parameters['Ki_comp1'].get()
        ki2 = self.parameters['Ki_comp2'].get()
        
        potencia_rel = ki2 / ki1 if ki1 > 0 else 0
        
        self.results_vars['ed50_comp1'].set(f"ED50 Comp1: {ki1:.3f} mmol/L")
        self.results_vars['ed50_comp2'].set(f"ED50 Comp2: {ki2:.3f} mmol/L")
        self.results_vars['potencia_rel'].set(f"Potencia relativa: {potencia_rel:.1f}x")
        self.results_vars['ki_ratio'].set(f"Ratio Ki2/Ki1: {potencia_rel:.1f}")
    
    def update_plots(self):
        """Actualizar todos los gráficos"""
        try:
            self.update_inhibition_plot()
            self.update_kinetic_plot()
            self.update_comparison_plot()
            self.update_results()
        except Exception as e:
            messagebox.showerror("Error", f"Error al actualizar gráficos: {str(e)}")
    
    def on_parameter_change(self):
        """Callback para cambios en parámetros"""
        # Actualización automática con un pequeño retraso para evitar demasiadas actualizaciones
        if hasattr(self, '_update_timer'):
            self.root.after_cancel(self._update_timer)
        self._update_timer = self.root.after(100, self.update_plots)
    
    def on_compound_change(self):
        """Callback para cambio de compuesto"""
        self.update_plots()
    
    def bind_events(self):
        """Vincular eventos"""
        # Actualizar gráficos cuando cambien los parámetros
        for param in self.parameters.values():
            if isinstance(param, (tk.DoubleVar, tk.IntVar)):
                param.trace("w", lambda *args: self.on_parameter_change())
    
    def reset_parameters(self):
        """Resetear parámetros a valores por defecto"""
        defaults = {
            'mu_max': 0.95,
            'Ki_comp1': 0.12,
            'Ki_comp2': 10.0,
            'K_max': 1e9,
            'ka': 0.1,
            'ke': 0.05,
            'initial_conc': 0.2,
            'time_points': 100
        }
        
        for key, value in defaults.items():
            if key in self.parameters:
                self.parameters[key].set(value)
        
        self.selected_compound.set("compound1")
        self.update_plots()
    
    def export_data(self):
        """Exportar datos actuales"""
        try:
            kinetic_data = self.calculate_kinetic_data()
            filename = f"kinetic_data_{self.selected_compound.get()}.csv"
            kinetic_data.to_csv(filename, index=False)
            messagebox.showinfo("Éxito", f"Datos exportados a {filename}")
        except Exception as e:
            messagebox.showerror("Error", f"Error al exportar datos: {str(e)}")
    
    def run(self):
        """Ejecutar la aplicación"""
        self.root.mainloop()

# Función principal
def main():
    """Función principal para ejecutar la aplicación"""
    print("Iniciando Modelo Cinético de Inhibición Enzimática - Versión Interactiva")
    print("="*70)
    
    app = InteractiveKineticModel()
    app.run()

if __name__ == "__main__":
    main()