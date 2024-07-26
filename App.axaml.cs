using Avalonia;
using Avalonia.Controls.ApplicationLifetimes;
using Avalonia.Markup.Xaml;
using EffectsDemoDRM.ViewModels;
using EffectsDemoDRM.Views;
using Avalonia.LinuxFramebuffer;
using System;

namespace EffectsDemoDRM
{
    public partial class App : Application
    {
        public override void Initialize()
        {
            AvaloniaXamlLoader.Load(this);
        }

        public override void OnFrameworkInitializationCompleted()
        {
            if (ApplicationLifetime is ISingleViewApplicationLifetime control)
            {
                control.MainView = new MainWindow {
                    DataContext = new MainWindowViewModel(),
                };
            }

            Console.WriteLine("Hello Torizon!");
            base.OnFrameworkInitializationCompleted();
        }
    }
}
