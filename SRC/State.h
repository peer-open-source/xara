#pragma once

enum class State : int
  {
    None           =     0, // save? push?

    Init           =  0400,
    InitVeloc      =  0100,
    InitAccel      =  0200,
    InitTotal      =  0700,

    Past           =   040,
    PastVeloc      =   010,
    PastAccel      =   020,
    PastTotal      =   070,

    Pres           =    04,
    PresVeloc      =    01,
    PresAccel      =    02,
    PresTotal      =    07,

//  allTotal       = Init_all | Past_all | Pres_all, // 0777

    set_uid_onVeloc = 04000,
    set_gid_onVeloc = 02000,
    Load            = 01000, // sticky_bit

//  state_mask      = all_all | set_uid_on_exe | set_gid_on_exe | Load      , // 07777

    state_not_known = 0xffff,

    add_state       = 0x1000,
    remove_state    = 0x2000,
    symlink_state   = 0x4000
  };

  inline constexpr State
  operator&(State x, State y)
  {
    return static_cast<State>
      (static_cast<int>(x) & static_cast<int>(y));
  }

  inline constexpr State
  operator|(State x, State y)
  {
    return static_cast<State>
      (static_cast<int>(x) | static_cast<int>(y));
  }

  inline constexpr State
  operator^(State x, State y)
  {
    return static_cast<State>
      (static_cast<int>(x) ^ static_cast<int>(y));
  }

  inline constexpr State
  operator~(State x)
  {
    return static_cast<State>(~static_cast<int>(x));
  }

  inline State &
  operator&=(State & x, State y)
  {
    x = x & y;
    return x;
  }

  inline State &
  operator|=(State & x, State y)
  {
    x = x | y;
    return x;
  }

  inline State &
  operator^=(State & x, State y)
  {
    x = x ^ y;
    return x;
  }
