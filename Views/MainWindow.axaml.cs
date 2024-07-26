using Avalonia.Controls;
using Avalonia.Rendering;

namespace EffectsDemoDRM.Views;

public partial class MainWindow : UserControl
{
    public MainWindow()
    {
        InitializeComponent();

        this.ButtonStartShader.Click += (sender, e) =>
        {
            if(this.ButtonStartShader?.Content?.ToString() == "PLAY") {
                this.Shader.CanDraw = true;
                this.ButtonStartShader.Content = "STOP";
            }

            else if(this.ButtonStartShader?.Content?.ToString() == "STOP") {
                this.Shader.CanDraw = false;
                this.ButtonStartShader.Content = "PLAY";
            }
        };
    }
}
