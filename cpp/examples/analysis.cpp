#include <sdft/sdft.h>

#include <matplotlibcpp17/pyplot.h>
#include <NumCpp.hpp>

using sdft::SDFT;
namespace py = pybind11;

const double pi = std::acos(-1.);

nc::NdArray<double> phase(const nc::NdArray<double>& t, const nc::NdArray<double>& f)
{
  nc::NdArray<double> cs(t.shape());

  cs[0] = 0;

  for (size_t i = 1; i < t.shape().size(); ++i)
  {
    const auto dt = t[i] - t[i - 1];
    const auto dp = 2 * pi * f[i] * dt;
    cs[i] = dp + cs[i - 1];
  }

  return cs;
}

nc::NdArray<double> phase(const nc::NdArray<double>& t, const double f)
{
  return phase(t, nc::zeros<double>(t.shape()) + f);
}

int main()
{
  py::scoped_interpreter guard{};


  // 1) generate input signal

  const auto sr = 44100;  // sample rate in Hz
  const auto n = 1 * sr;  // number of samples
  const auto m = 1000;    // number of dft bins

  const auto t = nc::arange<double>(n) / double(sr);  // timestamps in seconds

  const auto f = 1000.;  // single frequency
  // const auto f = nc::linspace<double>(0., 1000., n);  // linear chirp
  // const auto f = nc::sin(pi * t * t[-1]) * 1000. + 1000.;  // frequency wave

  const auto x = nc::sin(phase(t, f));  // sample vector of shape (n)


  // 2) estimate output dft

  SDFT<double> sdft(m);  // create sdft plan

  nc::NdArray<std::complex<double>> dft(n, m);  // dft matrix of shape (n, m)

  sdft.sdft(n, x.data(), dft.data());  // perform sdft analysis


  // 3) plot spectrogram

  const auto db = 20. * nc::log10(nc::abs(dft));

  auto plot = matplotlibcpp17::pyplot::import();

  auto data = nc::pybindInterface::nc2pybind(db.transpose());

  plot.imshow(Args(data), Kwargs(
    "extent"_a = py::make_tuple(0, n / sr, 0, sr / 2),
    "origin"_a = "lower",
    "aspect"_a = "auto",
    "cmap"_a = "inferno",
    "interpolation"_a = "nearest"));

  auto cbar = plot.colorbar();

  plot.xlabel(py::make_tuple("s"));
  plot.ylabel(py::make_tuple("Hz"));

  // TODO: cbar.set_label("dB");

  plot.ylim(py::make_tuple(0, 5000));
  plot.clim(py::make_tuple(-120, 0));

  plot.show();


  return 0;
}
